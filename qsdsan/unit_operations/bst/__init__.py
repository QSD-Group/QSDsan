# -*- coding: utf-8 -*-
"""BioSTEAM-inherited unit operations for QSDsan."""

from ._abstract import Mixer, Splitter, FakeSplitter, ReversedSplitter
from ._compressor import IsothermalCompressor
from ._distillation import (
    BinaryDistillation,
    ShortcutColumn,
    MESHDistillation,
    AdiabaticMultiStageVLEColumn,
)
from ._facilities import ProcessWaterCenter
from ._flash import Flash
from ._heat_exchanging import HeatExchangerNetwork, HXprocess, HXutility
from ._pumping import Pump
from ._tank import Tank, MixTank, StorageTank

# Raw BioSTEAM units that QSDsan does not subclass (no QSDsan-specific behavior to
# add) but surfaces here so users can build systems through ``qsdsan`` without
# ``import biosteam``. These are plain re-exports of the BioSTEAM classes; per the
# architecture rule we re-export rather than copy/subclass. ``test_public_api``
# asserts the identity so a BioSTEAM rename trips a test instead of a user.
from biosteam import (
    IsenthalpicValve,
    Stripper,
    MolecularSieve,
    BatchBioreactor,
    VacuumSystem,
    Boiler,
    BoilerTurbogenerator,
    ChilledWaterPackage,
    CoolingTower,
    SolidsCentrifuge,
    )

__all__ = (
    'Mixer',
    'Splitter',
    'FakeSplitter',
    'ReversedSplitter',
    'IsothermalCompressor',
    'BinaryDistillation',
    'ShortcutColumn',
    'MESHDistillation',
    'AdiabaticMultiStageVLEColumn',
    'ProcessWaterCenter',
    'Flash',
    'HeatExchangerNetwork',
    'HXprocess',
    'HXutility',
    'Pump',
    'Tank',
    'MixTank',
    'StorageTank',
    # Raw BioSTEAM re-exports (not subclassed)
    'IsenthalpicValve',
    'Stripper',
    'MolecularSieve',
    'BatchBioreactor',
    'VacuumSystem',
    'Boiler',
    'BoilerTurbogenerator',
    'ChilledWaterPackage',
    'CoolingTower',
    'SolidsCentrifuge',
)
