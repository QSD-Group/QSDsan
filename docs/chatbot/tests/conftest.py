import sys
from pathlib import Path

# Make modules in docs/chatbot importable as top-level (import config, chunking, ...)
CHATBOT_DIR = Path(__file__).resolve().parents[1]
if str(CHATBOT_DIR) not in sys.path:
    sys.path.insert(0, str(CHATBOT_DIR))
