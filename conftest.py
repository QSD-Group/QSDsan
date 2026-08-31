import pytest


@pytest.fixture(autouse=True)
def _reset_doctest_state(request):
    """Isolate doctests from each other's global state. Clears the LCA registries,
    resets the auto-ID ticket counters, AND resets the default flowsheet so each
    doctest starts clean.

    The flowsheet reset matters under ``pytest-xdist``: doctests that read the
    default flowsheet (e.g. ``qsdsan._lca.LCA``, via ``create_example_system`` +
    ``Flowsheet.flowsheet.default.unit.M1``) would otherwise pick up units/streams
    left behind by another test sharing the worker, giving wrong impacts."""
    if not isinstance(request.node, pytest.DoctestItem):
        yield
        return
    import qsdsan as qs
    # `qs.default()` resets to a fresh 'default' flowsheet (clearing units/streams/
    # systems and the LCA registries), the utilities (incl. PowerUtility.price) and
    # CEPCI, and the auto-ID ticket counters (native + LCA).
    qs.default()
    yield
