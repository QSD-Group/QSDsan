"""FastAPI endpoint. Loads the index once at startup and wires real clients into
engine.answer_question. Kept thin: all logic is tested in engine/retrieval/prompts.
"""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, field_validator

import config
import embeddings
import engine
import retrieval

app = FastAPI(title="QSDsan docs chatbot")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # internal-first; tighten before any public launch
    allow_methods=["POST", "GET"],
    allow_headers=["*"],
)

_STATE = {"records": None, "claude": None}


class AskRequest(BaseModel):
    question: str

    @field_validator("question")
    @classmethod
    def _not_blank(cls, v):
        if not v or not v.strip():
            raise ValueError("question must not be blank")
        return v.strip()


class Citation(BaseModel):
    n: int
    url: str
    title: str
    source: str


class AskResponse(BaseModel):
    answer: str
    citations: list[Citation]


def _records():
    if _STATE["records"] is None:
        _STATE["records"] = retrieval.load_index(config.INDEX_URL)
    return _STATE["records"]


def _claude():
    if _STATE["claude"] is None:
        import anthropic

        _STATE["claude"] = anthropic.Anthropic()  # reads ANTHROPIC_API_KEY
    return _STATE["claude"]


def _answer(question: str, **kwargs):
    """Indirection seam so tests can monkeypatch the whole pipeline."""
    return engine.answer_question(
        question,
        records=_records(),
        embed_fn=embeddings.embed_texts,
        claude_client=_claude(),
        top_k=config.TOP_K,
        threshold=config.SIMILARITY_THRESHOLD,
        gen_model=config.GEN_MODEL,
    )


@app.get("/health")
def health():
    return {"status": "ok"}


@app.post("/ask", response_model=AskResponse)
def ask(req: AskRequest):
    return _answer(req.question)
