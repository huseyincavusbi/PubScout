# PubScout 🔬

PubScout is an AI-powered agent that actively monitors PubMed for new scientific literature on a chosen topic. It goes beyond simple alerts by using a Retrieval-Augmented Generation (RAG) pipeline to provide deep, contextual analysis of how new papers relate to existing research, delivering insights directly to your Slack workspace.

## ✨ Core Features

**Automated Monitoring**: Periodically polls PubMed for new papers based on a specific search query.

**AI-Powered Analysis**: Uses a Large Language Model (LLM) to perform granular data extraction, identifying key findings, significance, study type, sample size, and more.

**Long-Term Memory**: Employs a vector database to build a persistent knowledge base, enabling the agent to understand how new research confirms, contradicts, or extends prior work.

**Interactive Q&A**: Allows you to "ask" natural language questions of your research library directly from Slack (e.g., `/pubscout_ask Based on the current papers, what strategies are being used to make mRNA vaccines more effective in solid tumors?`).

**Real-Time Slack Alerts**: Delivers well-formatted, structured summaries of new papers directly to a designated Slack channel.

**Robust & Resilient**: Features a Dead-Letter Queue (DLQ) to ensure no paper is ever lost due to processing errors.

## ⚙️ How It Works

PubScout operates with two primary, decoupled services:

### 1. The Analysis Pipeline (Automated)

This is the core, automated process that finds and analyzes new papers. It runs continuously in the background.

**Poller Script ➔ PubMed API ➔ New Paper Found ➔ Agent Server ➔ Vector DB (Retrieve Context) ➔ LLM (Analyze) ➔ Slack Alert + Update Vector DB**

### 2. The Interactive Query Pipeline

This workflow is triggered on-demand by a user in Slack, providing an interface to the knowledge base built by the analysis pipeline.

**User in Slack (`/pubscout_ask`) ➔ Interactive Bot ➔ Vector DB (Retrieve Context) ➔ LLM (Synthesize Answer) ➔ Post Reply in Slack**

## 🚀 Tech Stack

The project uses a modern stack chosen for performance, intelligence, and ease of development.

| Category | Technology | Why It Was Chosen |
|----------|------------|-------------------|
| AI Reasoning | Groq & Llama 4 Scout | Provides exceptionally fast and high-quality LLM inference, which is crucial for real-time analysis and interaction. |
| AI Memory | ChromaDB | A powerful and easy-to-use vector database that creates the project's long-term, searchable knowledge base. |
| AI Embeddings | Sentence-Transformers | A state-of-the-art library for converting text into meaningful numerical vectors, which powers the semantic search. |
| API Interaction | BioPython | The standard for interacting with the NCBI PubMed API, providing a reliable way to fetch publication data. |
| Notifications & UI | Slack Bolt & SDK | The Bolt framework makes it simple to build robust, interactive Slack commands, while the SDK is used for sending alerts. |
| Scheduling | Schedule | A lightweight and intuitive library for running the PubMed poller on a recurring schedule within the script. |

## 🛠️ Setup and Installation

Follow these steps to get PubScout running locally.

### 1. Clone the Repository

```bash
git clone https://github.com/huseyincavusbi/PubScout.git
cd PubScout
```

### 2. Create and Activate a Virtual Environment

```bash
# For Unix/macOS
python3 -m venv PubScout
source PubScout/bin/activate

# For Windows
python -m venv PubScout
.\PubScout\Scripts\activate
```

### 3. Install Dependencies

The project's dependencies are listed in requirements.txt.

```bash
pip install -r requirements.txt
```

### 4. Configure Your Secrets

An example environment file is provided. Create your own `.env` file from the example to store your secret API keys.

```bash
cp .env.example .env
```

Now, edit the `.env` file with your credentials:

```env
# PubScout Environment Variables
# See Slack API docs for tokens: https://api.slack.com/apps
# You need to enable Socket Mode for the interactive bot.

PUBMED_EMAIL="your.email@example.com"
GROQ_API_KEY="gsk_..."
SLACK_BOT_TOKEN="xoxb-..."
SLACK_APP_TOKEN="xapp-..."
SLACK_CHANNEL_ID="C12345678"
```

### 5. Configure the Agent

Open `config.py` and set the `SEARCH_QUERY` to the topic you want to track.

```python
# The topic you want to track. Be specific for good results.
SEARCH_QUERY = "mRNA Vaccines on Solid Tumors" # Or whatever your topic is
```

## ▶️ How to Run

PubScout consists of two primary services that should be run in separate terminals: the **Analysis Agent** (which processes new papers) and the **Interactive Scout** (which answers your questions).

*Prerequisite: Ensure you have completed all steps in the "Setup and Installation" section.*

### Terminal 1: Run the Main Analysis Agent

This service is the engine of the project. It runs the poller to find new papers, analyzes them, and populates the vector database.

1.  **Run the `poller.py` script first** to populate the queue with new papers.
    ```bash
    python poller.py
    ```
    **Expected Output:**
    ```
    INFO - PubScout Poller started.
    INFO - Performing initial search on startup...
    INFO - Polling PubMed for new papers with query: 'mRNA Vaccines on Solid Tumors'
    INFO - Found 20 new paper(s). Writing to queue file.
    INFO - Scheduling job to run every 60 minutes.
    ```

2.  **Now, start the main agent server** to process the queue.
    ```bash
    python agent_server.py
    ```
    **Expected Terminal Output:** The agent will start up, find the queue file, and begin processing the batch.
    ```
    INFO - All clients and models initialized successfully.
    INFO - PubScout Agent Server started.
    INFO - Found 20 new papers in queue. Creating lock file...
    INFO - Submitting 20 papers to the thread pool...
    INFO - --- Starting processing for PMID: 39938154 ---
    ... (processing continues for all papers) ...
    INFO - --- Batch Processing Summary ---
    INFO - Total papers in batch: 20
    INFO - Successfully processed: 20
    INFO - Failed and sent to DLQ: 0
    INFO - -----------------------------
    INFO - Batch finished. Cleaning up queue and lock files.
    ```

> **📢 Note on Slack Alerts:** In parallel with the terminal output, you will see alerts for each successfully processed paper appear in your configured Slack channel. They will look like this:

<div align="center">
  <img src="Example Images/sent_slack.png" alt="Slack Alert Example" width="100%">
  <p><em>Example of automated Slack alerts for new research papers</em></p>
</div>

---

### Terminal 2: Run the Interactive Scout Bot

This service connects to Slack and listens for your commands, allowing you to query the knowledge base built by the Analysis Agent.

1.  **Start the interactive scout.**
    ```bash
    python interactive_scout.py
    ```
    **Expected Terminal Output:**
    ```
    INFO - Interactive Scout initialized successfully.
    INFO - Starting PubScout Interactive Bot...
    INFO - ⚡️ Bolt app is running!
    ```

## 💬 Interact with Your Scout in Slack

With both services running, you can now query the agent. The responses will depend entirely on the context available in its knowledge base.

### ✅ Example of a Successful Answer
When the agent finds relevant context, it will synthesize a direct answer.

**Command:** `/pubscout_ask Based on the current papers, what strategies are being used to make mRNA vaccines more effective in solid tumors?`

<div align="center">
  <img src="Example Images/success.png" alt="Successful Query Example" width="100%">
  <p><em>Example of a comprehensive answer based on available research context</em></p>
</div>

---

### 🛡️ Example of a Grounded Response (No Hallucination)
If the agent does not find specific context to answer your question, it will state so clearly instead of making up an answer. This is a critical feature that prevents misinformation.

<div align="center">
  <img src="Example Images/failure.png" alt="Grounded Response Example" width="100%">
  <p><em>Example of honest response when relevant context is not available</em></p>
</div>

---

## 🎯 Key Benefits

- **🤖 AI-Powered Intelligence**: Leverages cutting-edge LLMs for deep research analysis
- **🔍 Contextual Understanding**: Builds knowledge over time to provide meaningful insights
- **⚡ Real-Time Monitoring**: Never miss important developments in your field
- **🚫 No Hallucination**: Provides honest, grounded responses based on actual research
- **💬 Natural Language Interface**: Query your research library using everyday language
- **🔄 Continuous Learning**: Knowledge base grows with each new paper processed

## 📝 License

This project is open source and available under the [Apache License 2.0](LICENSE).

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📧 Contact

For questions or support, please open an issue on GitHub.
