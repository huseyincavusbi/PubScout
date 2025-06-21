# config.py
# Central configuration file for the PubScout agent.
# All non-secret settings should be defined here.

# --- PubScout Search Configuration ---
# The topic you want to track. Be specific for good results.
SEARCH_QUERY = "mRNA Vaccines on Solid Tumors" 

# --- PATHS ---
DB_PATH = "db/processed_papers.db"
NEW_PAPERS_FILE = "data/new_papers.json"
DEAD_LETTER_QUEUE_FILE = "data/failed_pmids.json" 
LOCK_FILE = "data/processing.lock"
LOG_FILE_POLLER = "logs/poller.log"
LOG_FILE_AGENT = "logs/agent.log"

# --- VECTOR DATABASE SETTINGS ---
VECTOR_DB_PATH = "db/vector_db"
# This is the model used to create vector embeddings. It's small, fast, and effective.
EMBEDDING_MODEL = 'all-MiniLM-L6-v2'

# --- API & MODEL ---
# The model used for analysis via the Groq API.
GROQ_MODEL = "meta-llama/llama-4-scout-17b-16e-instruct"

# --- POLLER SETTINGS ---
# How often the poller checks PubMed, in minutes.
POLL_INTERVAL_MINUTES = 60

# --- CONCURRENCY & PERFORMANCE ---
# The maximum number of papers to process at the same time.
MAX_CONCURRENT_WORKERS = 5
# The delay in seconds between starting each new paper processing job.
# Set to 2.0 to respect the 30 RPM rate limit of the Groq API (60s / 2s = 30 RPM).
REQUEST_PACING_SECONDS = 2.0

# --- LLM PROMPT ENGINEERING ---.
# Instructions for the Groq LLM to perform granular data extraction.
GRANULAR_ANALYSIS_PROMPT = """
Analyze the following biomedical abstract. Your only output must be a single, raw, valid JSON object.
Do not include any other text, explanations, or markdown formatting like ```json.
The JSON object must contain these exact keys:
- "one_sentence_summary": A single, concise summary sentence.
- "key_findings": A list of strings detailing the main results.
- "significance": A paragraph explaining the importance of the findings.
- "relation_to_prior_work": Based on the provided context, state if this paper confirms, contradicts, extends, or is unrelated to prior findings.
- "study_type": The type of study (e.g., "Randomized Controlled Trial", "Meta-Analysis", "Case Study", "In Vitro", "Not specified").
- "sample_size": The primary sample size as a single integer ONLY. If a number is not mentioned, if it's a range, or if it's not applicable, the value MUST be 0. Do not include any text.
- "key_molecules": A list of the most important molecules, proteins, or pathways mentioned. If none, provide an empty list [].

CONTEXT FROM PREVIOUSLY ANALYZED PAPERS:
---
{context}
---

NEW PAPER TO ANALYZE:
Title: {title}
Abstract: {abstract}
"""

# Instructions for the Groq LLM to synthesize an answer from retrieved context.
SYNTHESIS_PROMPT = """
You are PubScout, a helpful AI research assistant. Your task is to answer the user's question based *only* on the provided context from previously analyzed papers.
Do not use any external knowledge. If the context does not contain the answer, state that you cannot answer based on the available information.
Be concise and clear in your answer.

CONTEXT FROM RELEVANT PAPERS:
---
{context}
---

USER'S QUESTION:
{question}

Based on the context provided, what is the answer to the user's question?
"""