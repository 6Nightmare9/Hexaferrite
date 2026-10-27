"""Shared test helpers (importable as `_helpers` from any test: a separate module, not `conftest`, because
pytest may have several conftest.py files (tests/theory/conftest.py) and `from conftest import ...` is then ambiguous)."""


def run_as(modname: str, src: str):
    """Execute source as if it were a module named `modname` (to test role checks); return its scope."""
    scope = {"__name__": modname}
    exec(compile(src, f"<{modname}>", "exec"), scope)
    return scope
