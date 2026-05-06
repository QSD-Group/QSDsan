# QSDsan Unit Structure Design

Date: 2026-05-06

## Goal

QSDsan should make its unit operation library easier to understand by separating units into three behavior-based groups:

1. BioSTEAM-compatible units
2. Static QSDsan units
3. Dynamic QSDsan units

The structure should preserve existing imports while giving users a clearer way to choose the right unit class.

## Proposed Public Structure

Expose three public namespaces under `qsdsan.sanunits`:

```text
qsdsan.sanunits
|-- bst
|-- static
`-- dynamic
```

Example usage:

```python
from qsdsan.sanunits import bst, static, dynamic

P1 = bst.Pump(...)
T1 = static.SepticTank(...)
R1 = dynamic.CSTR(...)
```

Legacy imports such as `qs.sanunits.Pump` should keep working during migration.

## Unit Groups

### BioSTEAM-Compatible Units

These units should behave like BioSTEAM units, with only the compatibility layer needed for use in QSDsan systems. They may handle QSDsan stream initialization, naming, and light integration with `SanUnit`, but they should not add sanitation-specific design assumptions.

Likely examples include `Mixer`, `Splitter`, `FakeSplitter`, `ReversedSplitter`, `Pump`, `Flash`, `BinaryDistillation`, `ShortcutColumn`, `MESHDistillation`, `AdiabaticMultiStageVLEColumn`, `HXprocess`, `HXutility`, `HeatExchangerNetwork`, `StorageTank`, `MixTank`, `IsothermalCompressor`, and `ProcessWaterCenter`.

### Static QSDsan Units

These units represent sanitation, wastewater, and resource recovery technologies using steady-state assumptions. They may include QSDsan-specific design, costing, construction, equipment, TEA, LCA, and `WasteStream` behavior, but they are not primarily dynamic models.

Likely examples include sanitation fixtures, treatment beds, tanks, septic systems, sedimentation, screening, sludge treatment, sludge pasteurization, crop application, hydrothermal units, hydroprocessing units, electrochemical cells, membrane gas extraction, membrane distillation, and system-specific public units.

### Dynamic QSDsan Units

These units are the core QSDsan differentiator. They should expose dynamic simulation contracts consistently:

- `isdynamic` support
- `state` and `dstate`
- `_init_state`
- `_compile_AE` or `_compile_ODE`
- consistent state headers
- process-model coupling where relevant
- tests that compare against expected static or benchmark behavior

Likely examples include dynamic influent, hydraulic delay, dynamic mixers and splitters where applicable, CSTR/PFR-style reactors, suspended growth bioreactors, ASM/ADM-coupled units, membrane bioreactors with dynamic process models, and junctions for dynamic process interfaces.

## Migration Strategy

1. Add `qsdsan.sanunits.bst`, `qsdsan.sanunits.static`, and `qsdsan.sanunits.dynamic` as re-export namespaces.
2. Classify the existing unit classes into the three namespaces without moving implementation files first.
3. Keep existing `qsdsan.sanunits.<ClassName>` imports as compatibility aliases.
4. Add documentation explaining how users should choose among the three groups.
5. Add or update tests to verify that each namespace exports the intended classes.
6. Only move source files after the re-export structure is stable and tested.

## Design Rules

- Classify by behavior, not source package.
- Use `bst` only for classes that preserve BioSTEAM semantics.
- Put steady-state QSDsan assumptions in `static`.
- Reserve `dynamic` for units with explicit dynamic state behavior.
- Avoid duplicating BioSTEAM implementations unless QSDsan adds meaningful behavior.
- Prefer compatibility wrappers over copied code for imported BioSTEAM units.
- Preserve existing public APIs until a deprecation policy is written.

## Open Decisions

- Whether `qsdsan.sanunits.bst` should contain all wrapped BioSTEAM units or only those already tested with QSDsan streams.
- Whether dynamically capable classes should appear in both `static` and `dynamic` if they can run either way, or only in `dynamic` with static operation documented as a mode.
- Whether to create a formal deprecation timeline for direct imports from `qsdsan.sanunits`.
