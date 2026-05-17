# TODO

## LCA refactoring follow-ups

- **`_transportation_specs`**: `SanUnit` gained `_construction_specs` for declaring default construction materials as class-level data (resolved lazily by `LCA`). A parallel `_transportation_specs` mechanism for default transportation activities has not been implemented yet but is likely needed for consistency. When added, it should follow the same pattern: a class-level tuple of dicts with keys `item`, `load`, `load_unit`, `distance`, `distance_unit`, `interval`, `interval_unit`, resolved by `LCA._resolve_transportation_specs()` at LCA creation time.
