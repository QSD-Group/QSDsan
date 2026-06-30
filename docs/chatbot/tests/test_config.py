import importlib


def test_defaults_present():
    config = importlib.import_module("config")
    assert config.GEN_MODEL == "claude-haiku-4-5-20251001"
    assert config.EMBED_MODEL == "voyage-4-lite"
    assert config.TOP_K == 6
    assert 0.0 < config.SIMILARITY_THRESHOLD < 1.0
    assert config.QSDSAN_DOCS_BASE.startswith("https://qsdsan.readthedocs.io")
    assert config.EXPOSAN_RAW_BASE.startswith("https://raw.githubusercontent.com/QSD-Group/EXPOsan")
    assert config.EXPOSAN_BLOB_BASE.startswith("https://github.com/QSD-Group/EXPOsan")


def test_env_override(monkeypatch):
    monkeypatch.setenv("CHATBOT_GEN_MODEL", "claude-sonnet-4-6")
    import config
    importlib.reload(config)
    assert config.GEN_MODEL == "claude-sonnet-4-6"
