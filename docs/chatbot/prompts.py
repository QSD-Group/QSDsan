"""System/user prompts and guardrail text. Behavior of the bot lives here."""
import re

SYSTEM_PROMPT = (
    "You are the QSDsan and EXPOsan documentation assistant. "
    "Answer ONLY using the numbered excerpts provided in the user message. "
    "Cite every claim with its excerpt number and URL, like [1]. "
    "If the excerpts do not contain the answer, say you could not find it in the "
    "QSDsan/EXPOsan docs and do not guess. "
    "Refuse questions unrelated to QSDsan or EXPOsan. "
    "When you include code, keep it minimal and faithful to the excerpts. "
    "Format answers in GitHub-flavored Markdown: use ## headings for sections, "
    "**bold** for key terms, hyphen bullet lists, and fenced ``` code blocks for code."
)

CODE_DISCLAIMER = "Draft from the docs - verify against the linked pages."


def build_user_prompt(question: str, retrieved: list[dict]) -> str:
    """Render numbered excerpts followed by the question."""
    lines = ["Excerpts:"]
    for i, chunk in enumerate(retrieved, start=1):
        lines.append(
            f"[{i}] ({chunk['source']}) {chunk['title']} - {chunk['url']}\n{chunk['text']}"
        )
    lines.append("")
    lines.append(f"Question: {question}")
    return "\n".join(lines)


def append_code_disclaimers(answer_md: str) -> str:
    """Append the verify-against-docs note right after each fenced code block."""
    pattern = re.compile(r"```.*?```", re.DOTALL)

    def repl(match):
        return match.group(0) + f"\n\n*{CODE_DISCLAIMER}*"

    return pattern.sub(repl, answer_md)


def refusal_message() -> str:
    """Out-of-scope / low-similarity refusal that also explains what the bot does.

    Leading with the capability summary means meta questions ("what can you help
    me with?") get a useful answer, while still signalling the strict grounding.
    """
    return (
        "I'm the QSDsan and EXPOsan documentation assistant. I answer questions "
        "grounded in the QSDsan API and tutorials and the EXPOsan example systems, "
        "such as how to create a WasteStream, run the BSM1 system, or set up a "
        "dynamic simulation. I couldn't find this in the QSDsan/EXPOsan docs, so "
        "try rephrasing or use the search function at the top of the navigation bar."
    )
