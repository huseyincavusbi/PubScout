# poller.py 
# This script periodically checks PubMed for new papers based on the search
# query in config.py. It uses tenacity for robust network calls and professional
# logging for monitoring.

import os
import json
import time
import schedule
import sqlite3
import logging
import config  # Import our new config file
from Bio import Entrez
from dotenv import load_dotenv
from tenacity import retry, stop_after_attempt, wait_fixed

# --- CONFIGURATION ---
load_dotenv()
Entrez.email = os.getenv("PUBMED_EMAIL")

# --- LOGGING SETUP ---
os.makedirs(os.path.dirname(config.LOG_FILE_POLLER), exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(config.LOG_FILE_POLLER),
        logging.StreamHandler()
    ]
)

# --- DATABASE LOGIC ---
def setup_database():
    """Creates the SQLite database and table if they don't exist."""
    os.makedirs(os.path.dirname(config.DB_PATH), exist_ok=True)
    with sqlite3.connect(config.DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute("CREATE TABLE IF NOT EXISTS processed_papers (pmid TEXT PRIMARY KEY)")
        conn.commit()

def get_processed_pmids():
    """Fetches all PMIDs that we have already processed from the database."""
    if not os.path.exists(config.DB_PATH):
        return set()
    with sqlite3.connect(config.DB_PATH) as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT pmid FROM processed_papers")
        return {row[0] for row in cursor.fetchall()}

# --- CORE POLLING LOGIC ---
@retry(stop=stop_after_attempt(3), wait=wait_fixed(10))
def search_pubmed():
    """
    Checks PubMed for new papers, avoiding runs if the agent is busy.
    This function will automatically retry on network errors.
    """
    if os.path.exists(config.LOCK_FILE):
        logging.warning("Agent is currently processing. Skipping this polling cycle.")
        return

    logging.info(f"Polling PubMed for new papers with query: '{config.SEARCH_QUERY}'")
    try:
        handle = Entrez.esearch(db="pubmed", term=config.SEARCH_QUERY, sort="date", retmax="20")
        record = Entrez.read(handle)
        handle.close()

        found_pmids = set(record["IdList"])
        processed_pmids = get_processed_pmids()
        new_pmids = list(found_pmids - processed_pmids)

        if not new_pmids:
            logging.info("No new papers found.")
            return

        logging.info(f"Found {len(new_pmids)} new paper(s). Writing to queue file.")
        
        os.makedirs(os.path.dirname(config.NEW_PAPERS_FILE), exist_ok=True)
        with open(config.NEW_PAPERS_FILE, "w") as f:
            json.dump({"new_pmids": new_pmids}, f)

    except Exception as e:
        logging.error("An error occurred during the PubMed search.", exc_info=True)
        raise  # Re-raise the exception to trigger the tenacity retry

# --- SCHEDULER ---
if __name__ == "__main__":
    setup_database()
    logging.info("PubScout Poller started.")
    
    logging.info("Performing initial search on startup...")
    try:
        search_pubmed()
    except Exception as e:
        logging.critical("Initial PubMed search failed after multiple retries. Please check network/config.")

    logging.info(f"Scheduling job to run every {config.POLL_INTERVAL_MINUTES} minutes.")
    schedule.every(config.POLL_INTERVAL_MINUTES).minutes.do(search_pubmed)

    while True:
        schedule.run_pending()
        time.sleep(1)