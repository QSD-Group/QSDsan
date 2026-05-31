"""Run the seeded eval questions against a live /ask endpoint.

Usage:
    python docs/chatbot/eval/run_eval.py --endpoint http://localhost:8000/ask

Exits non-zero if any case fails, so it can gate a release manually. Requires a
running server with real keys; not part of the unit test suite.
"""
import argparse
import sys
from pathlib import Path

import requests
import yaml

QUESTIONS = Path(__file__).parent / "questions.yaml"


def run(endpoint: str) -> int:
    items = yaml.safe_load(QUESTIONS.read_text())
    failures = 0
    for item in items:
        resp = requests.post(endpoint, json={"question": item["question"]}, timeout=120)
        data = resp.json()
        answer = data.get("answer", "")
        citations = data.get("citations", [])
        if item["expect"] == "refuse":
            ok = "couldn't find this" in answer and not citations
        else:
            urls = " ".join(c["url"] for c in citations)
            ok = item["value"] in urls
        status = "PASS" if ok else "FAIL"
        if not ok:
            failures += 1
        print(f"[{status}] {item['question']}")
    print(f"\n{len(items) - failures}/{len(items)} passed")
    return 1 if failures else 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--endpoint", default="http://localhost:8000/ask")
    args = parser.parse_args()
    sys.exit(run(args.endpoint))
