"""Central configuration for the QSDsan/EXPOsan docs chatbot.

Every tunable is read from an environment variable with a safe default so the
same code runs locally, on Render, and in the readthedocs build.
"""
import os

# --- Models ---
GEN_MODEL = os.getenv("CHATBOT_GEN_MODEL", "claude-haiku-4-5-20251001")
EMBED_MODEL = os.getenv("CHATBOT_EMBED_MODEL", "voyage-4-lite")

# --- Retrieval ---
TOP_K = int(os.getenv("CHATBOT_TOP_K", "6"))
SIMILARITY_THRESHOLD = float(os.getenv("CHATBOT_SIMILARITY_THRESHOLD", "0.4"))

# --- Source locations ---
QSDSAN_DOCS_BASE = os.getenv(
    "CHATBOT_QSDSAN_DOCS_BASE", "https://qsdsan.readthedocs.io/en/beta/"
)
EXPOSAN_RAW_BASE = os.getenv(
    "CHATBOT_EXPOSAN_RAW_BASE",
    "https://raw.githubusercontent.com/QSD-Group/EXPOsan/main/exposan/",
)
EXPOSAN_BLOB_BASE = os.getenv(
    "CHATBOT_EXPOSAN_BLOB_BASE",
    "https://github.com/QSD-Group/EXPOsan/blob/main/exposan/",
)

# Where the published index lives; the server loads this on startup.
INDEX_URL = os.getenv(
    "CHATBOT_INDEX_URL",
    "https://qsdsan.readthedocs.io/en/beta/_static/chatbot/index.json",
)

# Pointers used when a question is about EXPOsan / the catalog of example systems.
# (Targets the beta docs version; switch to /en/latest/ when promoted to main.)
SYSTEMS_DOC_URL = os.getenv(
    "CHATBOT_SYSTEMS_DOC_URL",
    "https://qsdsan.readthedocs.io/en/beta/systems/index.html",
)
EXPOSAN_REPO_URL = os.getenv(
    "CHATBOT_EXPOSAN_REPO_URL",
    "https://github.com/QSD-Group/EXPOsan",
)
