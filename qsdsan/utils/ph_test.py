#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 08:30:34 2023

@author: jiananfeng
"""

from qsdsan.utils.ph import pH_solver

pH_solver(kw=10**-14, activity=True, check_precipitation = True,
          weak_chemicals={'PO4':0.05, 'NH3':0.05, 'arsenic':0.05, 'HF':0.05, 'lactic':0.05, 'ascorbic': 0.05, 'acetate':0.05, 'propionic':0.05, 'SO4':0.05},
          other_chemicals={(('Cl', -1),): (0.05,), (('Mg', 2),): (0.05,), (('Fe', 3),): (0.05,)})

pH_solver(kw=10**-14, activity=True, check_precipitation = True,
          weak_chemicals={'PO4':0.05,},
          other_chemicals={(('Mg', 2),): (0.05,), (('Fe', 3),): (0.05,)})

chemicals={(('H3PO4', 0), ('H2PO4', -1), ('HPO4', -2), ('PO4', -3)): (0.00752, 6.230000000000001e-08, 4.8e-13, 0.05), (('NH4', 1), ('NH3', 0)): (5.6e-10, 0.05), (('H3AsO4', 0), ('H2AsO4', -1), ('HAsO4', -2), ('AsO4', -3)): (0.005, 8e-08, 6e-10, 0.05), (('HF', 0), ('F', -1)): (0.00072, 0.05), (('HCH3H5O3', 0), ('CH3H5O3', -1)): (0.000138, 0.05), (('H2C6H6O6', 0), ('HC6H6O6', -1), ('C6H6O6', -2)): (7.900000000000001e-05, 1.6e-12, 0.05), (('CH3COOH', 0), ('CH3COO', -1)): (1.76e-05, 0.05), (('CH3CH2COOH', 0), ('CH3CH2COO', -1)): (1.3400000000000002e-05, 0.05), (('HSO4', -1), ('SO4', -2)): (0.012, 0.05), (('Cl', -1),): (0.05,), (('Mg', 2),): (0.05,)}

chemicals={(('H2CO3', 0), ('HCO3', -1), ('CO3', -2), ('PO4', -3)): (4.3*10**-7, 4.8*10**-11, 0.001), (('Ca', 2),): (0.001,)}
