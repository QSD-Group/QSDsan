#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Cheung <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''
from chemicals.elements import molecular_weight

__all__ = ('cod_test_stoichiometry', 'electron_acceptor_cod')



#%%

dichromate_oxidizing_elements = ('C', 'H', 'O', 'N', 'S', 'P')
dichromate_oxidizing_elements_set = frozenset(dichromate_oxidizing_elements)

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
    nC, nH, nO, nN, nS, nP = 0., 0., 0., 0., 0., 0.
    ne = - charge
    
    if 'C' in atoms:
        nC = atoms['C']
    if 'H' in atoms:
        nH = atoms['H']
    if 'O' in atoms:
        nO = atoms['O']
    if 'N' in atoms:
        nN = atoms['N']
    if 'S' in atoms:
        nS = atoms['S']
    if 'P' in atoms:
        nP = atoms['P']
    
    if nC <= 0 or nH <= 0:
        if not (len(atoms) == 1 and nH == 2):
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