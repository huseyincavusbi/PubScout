# agent_server.py
# PubScout Agent Server
# This script monitors for new PubMed paper IDs, fetches their details,
# performs granular analysis using a LLM, sends summaries to Slack,
# stores embeddings in a vector database, and tracks processed papers.
# It uses a Dead-Letter Queue (DLQ) to ensure failed papers are not lost,
# and provides a summary report after each processing batch.

import os
import re
import json
import time
import sqlite3
import logging
import threading  # For creating a thread-safe lock for the DLQ
import config
from Bio import Entrez
from dotenv import load_dotenv
from groq import Groq, RateLimitError, InternalServerError
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError
from tenacity import retry, stop_after_attempt, wait_random_exponential, retry_if_exception_type
from concurrent.futures import ThreadPoolExecutor, as_completed 

# --- RAG IMPORTS ---
import chromadb
from sentence_transformers import SentenceTransformer

# --- CONFIGURATION & CLIENTS ---
load_dotenv()
GROQ_API_KEY = os.getenv("GROQ_API_KEY")
SLACK_BOT_TOKEN = os.getenv("SLACK_BOT_TOKEN")
SLACK_CHANNEL_ID = os.getenv("SLACK_CHANNEL_ID")
PUBMED_EMAIL = os.getenv("PUBMED_EMAIL")

# --- LOGGING SETUP ---
os.makedirs(os.path.dirname(config.LOG_FILE_AGENT), exist_ok=True)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(threadName)s - %(levelname)s - %(message)s', handlers=[logging.FileHandler(config.LOG_FILE_AGENT), logging.StreamHandler()])

# --- INITIALIZE CLIENTS & MODELS ---
try:
    groq_client = Groq(api_key=GROQ_API_KEY)
    slack_client = WebClient(token=SLACK_BOT_TOKEN)
    Entrez.email = PUBMED_EMAIL
    
    os.makedirs(os.path.dirname(config.VECTOR_DB_PATH), exist_ok=True)
    chroma_client = chromadb.PersistentClient(path=config.VECTOR_DB_PATH)
    embedding_model = SentenceTransformer(config.EMBEDDING_MODEL)
    vector_collection = chroma_client.get_or_create_collection(name="pubmed_papers")
    
    logging.info("All clients and models initialized successfully.")
except Exception:
    logging.critical("Failed to initialize clients or models.", exc_info=True)
    exit()

# --- GLOBAL LOCK FOR DLQ ---
dlq_lock = threading.Lock()

# --- HELPER FUNCTIONS ---
def sanitize_for_slack_mrkdwn(text: str) -> str:
    if not isinstance(text, str): text = str(text)
    return text.replace('&', '&').replace('<', '<').replace('>', '>')

# --- TOOL DEFINITIONS ---
def get_new_paper_ids():
    if not os.path.exists(config.NEW_PAPERS_FILE): return []
    try:
        with open(config.NEW_PAPERS_FILE, "r") as f: return json.load(f).get("new_pmids", [])
    except (FileNotFoundError, json.JSONDecodeError):
        logging.warning(f"Could not read or parse {config.NEW_PAPERS_FILE}.")
        return []

# --- DEAD-LETTER QUEUE FUNCTION ---
def add_to_dead_letter_queue(pmid: str):
    """Safely appends a failed PMID to the dead-letter queue JSON file."""
    with dlq_lock:
        try:
            os.makedirs(os.path.dirname(config.DEAD_LETTER_QUEUE_FILE), exist_ok=True)
            if os.path.exists(config.DEAD_LETTER_QUEUE_FILE):
                with open(config.DEAD_LETTER_QUEUE_FILE, "r") as f:
                    failed_data = json.load(f)
            else:
                failed_data = {"failed_pmids": []}
            
            if pmid not in failed_data["failed_pmids"]:
                failed_data["failed_pmids"].append(pmid)
            
            with open(config.DEAD_LETTER_QUEUE_FILE, "w") as f:
                json.dump(failed_data, f, indent=4)
            
            logging.info(f"PMID {pmid} has been added to the dead-letter queue.")
        except Exception:
            logging.critical(f"Could not write PMID {pmid} to the dead-letter queue!", exc_info=True)

@retry(wait=wait_random_exponential(min=2, max=30), stop=stop_after_attempt(3), retry=retry_if_exception_type(IOError))
def fetch_paper_details(pmid: str):
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml", retmode="text")
        records = Entrez.read(handle)
        handle.close()
        article = records['PubmedArticle'][0]['MedlineCitation']['Article']
        title = article.get('ArticleTitle', 'No Title Found')
        abstract = " ".join(article.get('Abstract', {}).get('AbstractText', []))
        return {"title": title, "abstract": abstract}
    except Exception:
        logging.error(f"Could not fetch details for PMID {pmid}. Retrying...", exc_info=True)
        raise

@retry(wait=wait_random_exponential(min=5, max=60), stop=stop_after_attempt(5), retry=retry_if_exception_type((RateLimitError, InternalServerError)))
def analyze_with_groq(title: str, abstract: str, context: str = "No prior context available."):
    prompt = config.GRANULAR_ANALYSIS_PROMPT.format(context=context, title=title, abstract=abstract)
    try:
        logging.info(f"Sending request to Groq API for granular analysis for '{title[:50]}...'")
        chat_completion = groq_client.chat.completions.create(messages=[{"role": "user", "content": prompt}], model=config.GROQ_MODEL, response_format={"type": "json_object"})
        return json.loads(chat_completion.choices[0].message.content)
    except (RateLimitError, InternalServerError) as e:
        logging.warning(f"Groq API error encountered: {type(e).__name__}. Tenacity will retry.")
        raise
    except Exception:
        logging.error("A non-retriable error occurred with the Groq API call.")
        raise

def send_slack_alert(pmid: str, title: str, analysis: dict):
    clean_title = re.sub(r'<.*?>', '', title)
    safe_header_title = (clean_title[:147] + "...") if len(clean_title) > 150 else clean_title
    study_type = analysis.get('study_type', 'N/A')
    sample_size = analysis.get('sample_size', 0)
    sample_size_text = str(sample_size) if sample_size > 0 else "N/A"
    molecules_list = analysis.get('key_molecules', [])
    molecules_text = ", ".join(molecules_list) if molecules_list else "N/A"
    blocks = [
        {"type": "header", "text": {"type": "plain_text", "text": safe_header_title}},
        {"type": "divider"},
        {"type": "section", "fields": [{"type": "mrkdwn", "text": f"*Study Type:*\n{study_type}"}, {"type": "mrkdwn", "text": f"*Sample Size:*\n{sample_size_text}"}, {"type": "mrkdwn", "text": f"*Key Molecules:*\n{molecules_text}"}]},
        {"type": "divider"},
        {"type": "section", "text": {"type": "mrkdwn", "text": f"*Summary:*\n{sanitize_for_slack_mrkdwn(analysis.get('one_sentence_summary', 'N/A'))}"}},
        {"type": "section", "text": {"type": "mrkdwn", "text": f"*Key Findings:*\n{'\n'.join([f'• {sanitize_for_slack_mrkdwn(item)}' for item in analysis.get('key_findings', [])]) if analysis.get('key_findings') else 'N/A'}"}},
        {"type": "section", "text": {"type": "mrkdwn", "text": f"*Significance:*\n{sanitize_for_slack_mrkdwn(analysis.get('significance', 'N/A'))}"}},
        {"type": "section", "text": {"type": "mrkdwn", "text": f"*Relation to Prior Work:*\n{sanitize_for_slack_mrkdwn(analysis.get('relation_to_prior_work', 'N/A'))}"}},
        {"type": "divider"},
        {"type": "context", "elements": [{"type": "mrkdwn", "text": f"Source: <https://pubmed.ncbi.nlm.nih.gov/{pmid}|PubMed ID: {pmid}>"}]}
    ]
    try:
        slack_client.chat_postMessage(channel=SLACK_CHANNEL_ID, blocks=blocks, text=f":scroll: New Publication: {clean_title}")
        logging.info(f"Successfully sent alert to Slack for PMID: {pmid}")
    except SlackApiError as e:
        logging.error(f"Failed to send Slack alert for PMID {pmid}. Reason: {e.response['error']}")

def mark_paper_as_processed(pmid: str):
    try:
        with sqlite3.connect(config.DB_PATH) as conn: conn.cursor().execute("INSERT OR IGNORE INTO processed_papers (pmid) VALUES (?)", (pmid,))
    except sqlite3.Error:
        logging.error(f"Could not write PMID {pmid} to database.", exc_info=True)

# --- WORKER FUNCTION WITH DLQ LOGIC ---
def process_single_paper(pmid: str) -> bool:
    """
    Defines the complete workflow for processing a single paper.
    This function is designed to be run in a separate thread.
    Returns True on success, False on failure.
    """
    logging.info(f"--- Starting processing for PMID: {pmid} ---")
    
    with sqlite3.connect(config.DB_PATH) as conn:
        if conn.cursor().execute("SELECT pmid FROM processed_papers WHERE pmid=?", (pmid,)).fetchone():
            logging.info(f"PMID {pmid} already processed. Skipping.")
            return True

    try:
        details = fetch_paper_details(pmid)
        if not details or not details.get("abstract"):
            logging.warning(f"No abstract for PMID {pmid}. Marking as processed and skipping.")
            mark_paper_as_processed(pmid)
            return True

        text_to_embed = f"Title: {details['title']}\nAbstract: {details['abstract']}"
        new_embedding = embedding_model.encode(text_to_embed).tolist()

        logging.info(f"Querying vector DB for papers similar to PMID {pmid}...")
        similar_results = vector_collection.query(query_embeddings=[new_embedding], n_results=3)
        
        context_for_prompt = "No similar papers found in the database yet."
        if similar_results and similar_results['ids'][0]:
            context_lines = []
            for i, pmid_sim in enumerate(similar_results['ids'][0]):
                meta = similar_results['metadatas'][0][i]
                context_lines.append(f"- PMID {pmid_sim} (Title: {meta.get('title', 'N/A')}): {meta.get('summary', 'Summary not available.')}")
            context_for_prompt = "\n".join(context_lines)
        logging.info(f"Found context for RAG prompt for PMID {pmid}.")

        analysis = analyze_with_groq(details["title"], details["abstract"], context=context_for_prompt)
        logging.info(f"Granular analysis successful for PMID: {pmid}")

        send_slack_alert(pmid, details["title"], analysis)

        logging.info(f"Adding PMID {pmid} to vector DB...")
        vector_collection.add(
            embeddings=[new_embedding],
            documents=[text_to_embed],
            metadatas=[{"pmid": pmid, "title": details['title'], "summary": analysis.get("one_sentence_summary")}],
            ids=[pmid]
        )
        
        mark_paper_as_processed(pmid)
        logging.info(f"--- Finished processing for PMID: {pmid} ---")
        return True 
    
    except Exception:
        # Error handling
        logging.error(f"A critical error occurred while processing PMID {pmid}. Adding to dead-letter queue.", exc_info=True)
        add_to_dead_letter_queue(pmid) # Call the DLQ function
        return False 

# --- MAIN AGENT LOOP ---
def main():
    """The main orchestration loop, now with failure tracking."""
    logging.info("PubScout Agent Server started.")
    
    while True:
        pmids = get_new_paper_ids()
        if not pmids:
            time.sleep(30)
            continue

        logging.info(f"Found {len(pmids)} new papers in queue. Creating lock file...")
        try:
            with open(config.LOCK_FILE, "w") as f: f.write(str(time.time()))
        except IOError:
            logging.critical("Could not create lock file.", exc_info=True)
            time.sleep(60)
            continue

        # --- CONCURRENT EXECUTION BLOCK ---
        futures = []
        with ThreadPoolExecutor(max_workers=config.MAX_CONCURRENT_WORKERS) as executor:
            logging.info(f"Submitting {len(pmids)} papers to the thread pool...")
            for pmid in pmids:
                future = executor.submit(process_single_paper, pmid)
                futures.append(future)
                time.sleep(config.REQUEST_PACING_SECONDS)
        
        # --- BATCH COMPLETION & ANALYSIS ---
        successful_count = 0
        failed_count = 0
        for future in as_completed(futures):
            if future.result():
                successful_count += 1
            else:
                failed_count += 1
        
        logging.info("--- Batch Processing Summary ---")
        logging.info(f"Total papers in batch: {len(pmids)}")
        logging.info(f"Successfully processed: {successful_count}")
        logging.info(f"Failed and sent to DLQ: {failed_count}")
        logging.info("-----------------------------")
        
        try:
            logging.info("Batch finished. Cleaning up queue and lock files.")
            if os.path.exists(config.NEW_PAPERS_FILE): os.remove(config.NEW_PAPERS_FILE)
            if os.path.exists(config.LOCK_FILE): os.remove(config.LOCK_FILE)
        except OSError:
            logging.error("Could not clean up files.", exc_info=True)

if __name__ == "__main__":
    main()