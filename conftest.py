import pytest


@pytest.fixture(autouse=True)
def _reset_lca_state(request):
    """Isolate LCA doctests: clear per-flowsheet LCA registries and reset the
    auto-ID ticket counters so each doctest starts from a clean, deterministic state."""
    if not isinstance(request.node, pytest.DoctestItem):
        yield
        return
    from qsdsan import ImpactIndicator, ImpactItem, Construction, Transportation
    for cls in (ImpactIndicator, ImpactItem, Construction, Transportation):
        cls.clear_registry(print_msg=False)
        cls.ticket_numbers.clear()
    yield
