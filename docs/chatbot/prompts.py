"""System/user prompts and guardrail text. Behavior of the bot lives here."""
import re

import config

SYSTEM_PROMPT = (
    "You are the QSDsan documentation assistant. "
    "Answer ONLY using the numbered excerpts provided in the user message. "
    "Cite every claim with its excerpt number and URL, like [1]. "
    "If the excerpts do not actually answer the question, do NOT guess and do NOT "
    "list generic suggestions. Reply in one or two sentences that the QSDsan "
    "documentation does not appear to cover it, and point the user to the search "
    "bar at the top of the page or to opening an issue on the QSDsan GitHub "
    "repository (https://github.com/QSD-Group/QSDsan). "
    "Refuse questions unrelated to QSDsan or EXPOsan. "
    "Excerpts are tagged by source: prefer (qsdsan) prose excerpts to EXPLAIN "
    "concepts, and prefer (code) excerpts (doctests and signatures from the "
    "source) for the EXACT, runnable usage and API. "
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
    """Technical question the docs do not cover: point to search + the GitHub repo."""
    return (
        "I couldn't find this in the QSDsan documentation. Try the search bar at the "
        "top of the page, or open an issue on the "
        f"[QSDsan GitHub repo]({config.QSDSAN_REPO_URL})."
    )


# Greetings and capability/meta questions. Matched against the whole (normalized)
# message, so a real question that merely contains one of these words is unaffected.
_SMALLTALK = {
    "hi", "hello", "hey", "yo", "hiya", "howdy", "hi there", "hello there",
    "good morning", "good afternoon", "good evening", "good day",
    "thanks", "thank you", "ty", "thx", "cheers",
    "how are you", "what's up", "whats up", "sup",
    "help", "what can you do", "what can you help with",
    "what can you help me with", "what do you do", "who are you",
    "what are you", "what is this", "what can i ask", "what can i ask you",
}


def is_smalltalk(question: str) -> bool:
    """True for greetings and capability/meta questions (non-technical)."""
    norm = question.strip().lower().strip(" !.?")
    return norm in _SMALLTALK


def greeting_message() -> str:
    """Friendly welcome for greetings and capability/meta questions."""
    return (
        "Hi! I'm the QSDsan documentation assistant. Ask me about the QSDsan API or "
        "tutorials, such as how to create a WasteStream, set up a dynamic simulation, "
        "or build a custom unit operation. What would you like to know?"
    )


# Catalog / discovery phrases: "what systems exist", "what have you built". These
# always get the pointer (the index is not a catalog of systems).
_EXPOSAN_CATALOG = (
    "what system",
    "which system",
    "list of system",
    "available system",
    "example system",
    "what have you built",
    "what did you build",
    "what can you simulate",
    "what models are available",
)

# Specific EXPOsan systems / the repo. Generic model names that also exist inside
# QSDsan (adm, asm, cas, hap, pm2) are deliberately excluded to avoid hijacking
# QSDsan questions.
_EXPOSAN_NAMES = (
    "exposan",
    "bsm1",
    "bsm2",
    "bwaise",
    "biobinder",
    "biogenic refinery",
    "new generator",
    "reclaimer",
    "eco_san",
    "eco-san",
    "eco san",
    "pou disinfection",
)

# Build/run/code intent. A specific system asked about this way is answerable from
# the indexed create_system wiring, so it should NOT be diverted to the pointer.
_CODE_INTENT = (
    "how",  # also catches "show me ..."
    "run",
    "load",
    "create",
    "build",
    "set up",
    "setup",
    "simulate",
    "import",
    "code",
    "example",
    "assemble",
    "construct",
)


def is_exposan_question(question: str) -> bool:
    """True when a question should be routed to the EXPOsan pointer response.

    Catalog/discovery questions always point to the Systems page + EXPOsan repo.
    A specific system asked about in a build/run/code way instead falls through
    to retrieval, so it can be answered from the indexed create_system wiring.
    """
    q = question.lower()
    if any(keyword in q for keyword in _EXPOSAN_CATALOG):
        return True
    if any(name in q for name in _EXPOSAN_NAMES):
        return not any(intent in q for intent in _CODE_INTENT)
    return False


def exposan_pointer_message() -> str:
    """Markdown pointer to the Systems page and the EXPOsan repository."""
    return (
        "For EXPOsan and the example systems that QSDsan has been used to build, "
        "these are the best places to look:\n\n"
        f"- QSDsan Systems page: [browse the systems]({config.SYSTEMS_DOC_URL})\n"
        f"- EXPOsan on GitHub: [QSD-Group/EXPOsan]({config.EXPOSAN_REPO_URL})"
    )
