#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
from chemicals.elements import molecular_weight
from thermosteam.reaction import (
    Reaction as Rxn,
    ParallelReaction as PRxn
    )
from . import auom

__all__ = (
    'cod_test_stoichiometry',
    'compute_stream_COD',
    'electron_acceptor_cod',
    'get_bmp_stoichiometry',
    'get_cod_stoichiometry',
    'get_digestion_rxns',
    )


#%%

dichromate_oxidizing_elements = ('C', 'H', 'O', 'N', 'S', 'P')
dichromate_oxidizing_elements_set = frozenset(dichromate_oxidizing_elements)

get_CHONSP = lambda atoms: \
    [atoms.get(atom) or 0. for atom in dichromate_oxidizing_elements]

def cod_test_stoichiometry(atoms, charge=0, MW=None, missing_handling='elemental'):
    r'''
    Return a dictionary of stoichiometric coefficients of the oxidation reaction
    by dichromate, given a dictionary of a molecule's constituent atoms and their
    counts, as well as the number of negative charge, if any.

    This function is based on the oxidation of organic materials by dichromate in
    an acid solution, as in a typical COD test of water or wastewater samples.
    Only C, H, O, N, S, P are considered active in the reaction.

    Parameters
    ----------
    atoms : dict[str, int or float]
        Dictionary of atoms and their counts, [-].
    charge : int or float
        Charge of the molecule.
    MW : float, optional
        Molecular weight of chemical, used only if `missing_handling` is 'Ash',
        [g/mol]
    missing_handling : str, optional
        How to handle compounds which do not appear in the stoichiometric
        reaction below. If 'elemental', return those atoms in the monatomic
        state; if 'ash', converts all missing attoms to 'Ash' in the output at
        a `MW` of 1 g/mol, [-]

    Returns
    -------
    stoichiometry : dict[str, float]
        Stoichiometric coefficients of the redox reaction. May inlcude the following
        keys for complete oxidation: 'H2O', 'CO2', 'NH4+', 'SO4-2', 'PO4-3';
        if `missing_handling` is 'elemental' can include the other elements;
        if `missing_handling` is 'ash', Ash will be present in the output
        if the compounds whose reactions are not included here.
        'Cr2O7-2' is always present, with negative values indicating dichromate is
        required/consumed. [-]

    .. note::

        The stoichiometry is given by:

            .. math::
                C_n H_a O_b N_c S_d P_e^{f-} + xCr_2O_7^{2-} + (8x+c-2d-3e+f)H^{+}
                    -> nCO_2 + 2xCr^{3+} + cNH_4^{+} + dSO_4^{2-} + ePO_4^{3-} + (b+7x-2n-4d-4e)H_2O
            .. math::
                x = \frac{4n+a-2b-3c+6d+5e+f}{6}

        Also included in the results is the moles of Cr2O7-2 required per mole of
        the mixture of the molecule.

        All products are in aqueous solution.

        Atoms not in ['C', 'H', 'N', 'O', 'S', 'P'] are returned as pure species;
        i.e. sodium hydroxide produces water and pure Na.

    Examples
    --------
    >>> # Acetate in COD test:
    >>> cod_test_stoichiometry({'C': 2, 'H':3, 'O':2}, -1)
    {'Cr2O7-2': -1.3333333333333333,
     'H+': -11.666666666666666,
     'Cr+3': 2.6666666666666665,
     'CO2': 2,
     'H2O': 7.333333333333332}
    '''
    products = {}
    nC, nH, nO, nN, nS, nP = get_CHONSP(atoms)
    ne = - charge

    if nC <= 0 or nH <= 0:
        if not (len(atoms) == 1 and nH == 2): # H2
            return {'Cr2O7-2': 0.}

    nCO2 = nC
    nNH4 = nN
    nSO4 = nS
    nPO4 = nP
    nCr2O7 = -(4*nC + nH - 2*nO - 3*nN + 6*nS + 5*nP + ne)/6
    nCr = -2*nCr2O7
    nH2O = nO - 7*nCr2O7 - 2*nC - 4*nS - 4*nP
    n_proton = 8*nCr2O7 - nN + 2*nS + 3*nP - ne

    if nCr2O7 != 0.0:
        products['Cr2O7-2'] = nCr2O7
    if n_proton != 0.0:
        products['H+'] = n_proton
    if nCr != 0.0:
        products['Cr+3'] = nCr
    if nCO2 !=  0.0:
        products['CO2'] = nCO2
    if nSO4 != 0.0:
        products['SO4-2'] = nSO4
    if nNH4 != 0.0:
        products['NH4+'] = nNH4
    if nPO4 != 0.0:
        products['PO4-3'] = nPO4
    if nH2O != 0.0:
        products['H2O'] = nH2O

    missing_handling = missing_handling.lower()
    if missing_handling == 'elemental':
        for atom, value in atoms.items():
            if atom not in dichromate_oxidizing_elements_set:
                products[atom] = value
    elif missing_handling == 'ash':
        cod_atoms = {i: atoms.get(i, 0) for i in dichromate_oxidizing_elements}
        MW = MW or molecular_weight(atoms)
        Ash = MW - molecular_weight(cod_atoms)
        if Ash/MW > 0.0001:
            products['Ash'] = Ash
    else:
        raise ValueError("Allowed values for `missing_handling` are 'elemental' and 'ash'.")
    return products


def electron_acceptor_cod(atoms, charge=0):
    if atoms == {'O':2}:
        return -1
    elif atoms == {'N':2}:
        return -1.5
    elif atoms == {'N':1, 'O':2} and charge == -1:
        return -1.5
    elif atoms == {'N':1, 'O':3} and charge == -1:
        return -2


def get_cod_stoichiometry(component, aqueous=False, **replace):
    r'''
    Get the molar stoichiometry for the theoretical
    chemical oxygen demand (COD) of a given component.

    COD stoichiometry is consistent with :func:`qsdsan.utils.cod_test_stoichiometry`
    other than the oxidant is O2 rather than Cr2O7-2,


    For organic components, elements other than "C", "H", "O", "N", "S", and "P" will
    be turned into "Ash" with a molecular weight of 1 g/mol.

    For inorganic components, all dict values will be 0.

    If `aqueous` == False, the stoichiometry is given by:

    .. math::
        C_nH_aO_bN_cS_dP_e + \frac{2n+0.5a-b-1.5c+3d+2.5e}{2}O_2
        -> nCO_2 + \frac{a-3c-2d}{2}H_2O + cNH_3 + dH_2SO_4 + \frac{e}{4}P_4O_{10}

    otherwise:

    .. math::
        C_nH_aO_bN_cS_dP_e + \frac{2n+0.5a-b-1.5c+3d+2.5e}{2}O_2 + (c-2d-3e)H^+
        -> nCO_2 + \frac{a-3c-2d-3e}{2}H_2O + cNH_4^+ + dSO_4^{2-} + ePO_4^{3-}

    Parameters
    ----------
    component : obj
        The component whose COD will be calculated.
    aqueous : bool
        Whether the reaction will happen in aqueous phase.
    replace : dict
        Alternative IDs of the reactant/product components,
        e.g., if S_O2 is the ID of dissolved oxygen instead of O2,
        then can pass replace={'O2': 'S_O2'}.

    Examples
    --------
    >>> from qsdsan import Component
    >>> from qsdsan.utils import get_cod_stoichiometry
    >>> Glucose = Component('Glucose', organic=True, particle_size='Soluble',
    ...                     degradability='Readily')
    >>> get_cod_stoichiometry(Glucose)
    {'Glucose': -1.0,
     'O2': -6.0,
     'CO2': 6,
     'H2O': 6.0,
     'NH3': 0.0,
     'H2SO4': 0.0,
     'P4O10': 0.0}
    '''
    cmp_ID = component.ID
    atoms = component.atoms

    keys = (cmp_ID, 'O2', 'CO2', 'H2O', 'NH3', 'H2SO4', 'P4O10') if not aqueous \
        else (cmp_ID, 'O2', 'H+', 'CO2', 'H2O', 'NH4+', 'SO42-', 'PO43-')
    dct = dict.fromkeys(keys, 0.)

    if atoms and component.organic:
        nC, nH, nO, nN, nS, nP = get_CHONSP(atoms)

        dct[cmp_ID] = -1.
        dct['O2'] = -(nC+0.25*nH-0.5*nO-0.75*nN+1.5*nS+1.25*nP)
        dct['CO2'] = nC
        if not aqueous:
            dct['H2O'] = 0.5*nH-1.5*nN-nS
            dct['NH3'] = nN
            dct['H2SO4'] = nS
            dct['P4O10'] = 0.25*nP
        else:
            dct['H+'] = -(nN-2*nS-3*nP)
            dct['H2O'] = 0.5*nH-1.5*nN-nS-1.5*nP
            dct['NH4+'] = nN
            dct['SO42-'] = nS
            dct['PO43-'] = nP

        cod_atoms = {i: atoms.get(i, 0) for i in dichromate_oxidizing_elements}
        MW = component.MW or molecular_weight(atoms)
        Ash = MW - molecular_weight(cod_atoms)
        if Ash/MW > 0.0001:
            dct['Ash'] = Ash

    for old_ID, new_ID in replace.items():
        dct[new_ID] = dct.pop(old_ID)

    return dct


def compute_stream_COD(stream, units='mg/L'):
    r'''
    Compute the chemical oxygen demand (COD) of a given stream
    by summing the COD of each component in the stream using:

    .. math::
        COD [\frac{kg}{m^3}] = mol_{component} [\frac{kmol}{m^3}] * \frac{g O_2}{mol\ component}
    '''
    try:
        COD = stream.COD
    except AttributeError: # not a WasteStream
        cmps = stream.components
        mol = stream.mol
        iCOD = np.array([-get_cod_stoichiometry(i)['O2'] for i in cmps])
        COD = (mol*iCOD).sum()*molecular_weight({'O': 2}) / stream.F_vol
    return auom('mg/L').convert(COD, units)


def get_bmp_stoichiometry(component, **replace):
    r'''
    Compute the theoretical biochemical methane potential (BMP) in
    mol :math:`CH_4`/mol chemical of a given component as in:

    .. math::
        C_vH_wO_xN_yS_z + \frac{4v-w-2x+3y+2z}{2}H2O ->
        \frac{4v+w-2x-3y-2z}{8}CH4 + \frac{(4v-w+2x+3y+2z)}{8}CO2 + yNH_3 + zH_2S

    For organic components, elements other than "C", "H", "O", "N", "S", and "P" will
    be turned into "Ash" with a molecular weight of 1 g/mol.

    For inorganic components, all dict values will be 0.

    Parameters
    ----------
    component : obj
        The component whose COD will be calculated.
    replace : dict
        Alternative IDs of the reactant/product components,
        e.g., if S_O2 is the ID of dissolved oxygen instead of O2,
        then can pass replace={'O2': 'S_O2'}.

    Examples
    --------
    >>> from qsdsan import Component
    >>> from qsdsan.utils import get_bmp_stoichiometry
    >>> Glucose = Component('Glucose', organic=True, particle_size='Soluble',
    ...                     degradability='Readily')
    >>> get_bmp_stoichiometry(Glucose)
    {'Glucose': -1.0, 'H2O': -0.0, 'CH4': 3.0, 'CO2': 3.0, 'NH3': 0.0, 'H2S': 0.0}
    '''
    cmp_ID = component.ID
    atoms = component.atoms
    keys = (cmp_ID, 'H2O', 'CH4', 'CO2', 'NH3', 'H2S')
    dct = dict.fromkeys(keys, 0.)

    if atoms and component.organic and component.formula != 'CH4':
        nC, nH, nO, nN, nS, nP = get_CHONSP(atoms)
        dct[cmp_ID] = -1.
        dct['H2O'] = -(nC-0.25*nH-0.5*nO+0.75*nN+0.5*nS)
        dct['CH4'] = 0.5*nC+0.125*nH-0.25*nO-0.375*nN-0.25*nS
        dct['CO2'] = 0.5*nC-0.125*nH+0.25*nO+0.375*nN+0.25*nS
        dct['NH3'] = nN
        dct['H2S'] = nS
        bmp_atoms = {i: atoms.get(i, 0) for i in dichromate_oxidizing_elements}
        MW = component.MW or molecular_weight(atoms)
        Ash = MW - molecular_weight(bmp_atoms)
        if Ash/MW > 0.0001:
            dct['Ash'] = Ash

    for old_ID, new_ID in replace.items():
        dct[new_ID] = dct.pop(old_ID)

    return dct


def get_digestion_rxns(components, X_biogas, X_growth, biomass_ID, biodegradability=1.):
    '''
    Generate anaerobic digestion (AD) and biomass growth reactions
    for a given set of components.

    AD stoichiometry is based on :func:`qsdsan.utils.get_bmp_stoichiometry`
    and biodegradabilities of the components as indicated in `biodegradability`.

    Biomass growth is purely based on mass balance, thus can potentially result
    in loss of atom balance.

    No reactions will be generated for inorganic components.

    Parameters
    ----------
    components : Iterable(obj)
        Set of components.
    X_biogas : float
        Fraction of the organic components that is used for AD.
    X_growth : float
        Fraction of the organic components that is used for biomass growth.
    biomass_ID : str
        ID of the biomass (should be included in the `components`).
    biodegradability : float or dict
        Biodegradabilities of the components.
        When given as a float, all organic components will be assumed to have the
        same biodegradability;
        when given as a dict, the keys should be the IDs of components and
        values the corresponding biodegradabilities,
        components without corresponding biodegradabilities will be assumed unbiodegradable.

    Examples
    --------
    >>> from qsdsan import Component, Components, set_thermo
    >>> from qsdsan.utils import load_example_cmps, get_digestion_rxns
    >>> example_cmps = load_example_cmps()
    >>> NH3 = Component('NH3', phase='g', organic=False, particle_size='Dissolved gas',
    ...                 degradability='Undegradable')
    >>> H2S = Component('H2S', phase='g', organic=False, particle_size='Dissolved gas',
    ...                 degradability='Undegradable')
    >>> P4O10 = Component('P4O10', phase='s',
    ...                   organic=False, particle_size='Particulate',
    ...                   degradability='Undegradable')
    >>> Biomass = Component('Biomass', phase='s', formula='CH1.8O0.5N0.2',
    ...                     organic=True, particle_size='Particulate',
    ...                     degradability='Slowly')
    >>> Ash = Component('Ash', phase='s', MW=1,
    ...                 organic=False, particle_size='Particulate',
    ...                 degradability='Undegradable')
    >>> for i in (P4O10, Biomass, Ash):
    ...    i.copy_models_from(example_cmps.NaCl, ('V',))
    >>> for i in (Biomass, Ash): i.default() # doctest:+ELLIPSIS
    {...
    >>> cmps = Components([*example_cmps, NH3, H2S, P4O10, Biomass, Ash])
    >>> cmps.compile()
    >>> set_thermo(cmps)
    >>> cmps
    CompiledComponents([H2O, CO2, N2O, NaCl, H2SO4, CH4, Methanol, Ethanol, NH3, H2S, P4O10, Biomass, Ash])
    >>> rxns = get_digestion_rxns(cmps, X_biogas=0.9, X_growth=0.07,
    ...                           biomass_ID='Biomass', biodegradability=0.87)
    >>> rxns # doctest: +SKIP
    ParallelReaction (by mol):
    index  stoichiometry                              reactant    X[%]
    [0]    Methanol -> 0.5 H2O + 0.25 CO2 + 0.75 CH4  Methanol   78.30
    [1]    Ethanol -> 0.5 CO2 + 1.5 CH4               Ethanol    78.30
    [2]    Methanol -> 1.3 Biomass                    Methanol    6.09
    [3]    Ethanol -> 1.87 Biomass                    Ethanol     6.09
    '''
    biomass = getattr(components, biomass_ID)
    biomass_MW = biomass.MW or molecular_weight(biomass.atoms)
    BD = dict.fromkeys(components.IDs, biodegradability) if isinstance(biodegradability, float) \
        else biodegradability

    if X_biogas+X_growth > 1:
        raise ValueError('Sum of `X_biogas`/`X_decomp` and `X_biogas` is '
                          f'{X_biogas+X_growth}, larger than 100%.')

    biogas_rxns, growth_rxns = [], []
    for i in components:
        ID = i.ID

        if ID == biomass_ID:
            continue

        X = BD.get(i.ID)
        if not X:
            continue # assume no entry means not biodegradable

        biogas_stoyk = get_bmp_stoichiometry(i)
        if not biogas_stoyk.get(i.ID): # no conversion of this chemical
            continue

        iX_biogas = X * X_biogas # the amount of component used for biogas production
        iX_growth = X * X_growth # the amount of component used for cell growth

        if iX_biogas:
            biogas_rxn = Rxn(reaction=biogas_stoyk, reactant=ID, X=iX_biogas,
                             check_atomic_balance=True)
            biogas_rxns.append(biogas_rxn)

        if iX_growth:
            growth_rxn = Rxn(f'{i.ID} -> {i.MW/biomass_MW}{biomass_ID}',
                             reactant=i.ID, X=iX_growth,
                             check_atomic_balance=False)
            growth_rxns.append(growth_rxn)

    if len(biogas_rxns)+len(growth_rxns)>1:
        return PRxn(biogas_rxns+growth_rxns)

    return []