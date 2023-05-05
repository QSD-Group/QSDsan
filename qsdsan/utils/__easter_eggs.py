# -*- coding: utf-8 -*-
"""
Created on Fri May  5 12:53:47 2023

@author: Yalin Li
"""

import webbrowser

prefix = 'https://qsdsan.readthedocs.io/en/latest/core_developers/'
suffix = '.html'

Joy = \
    'Joy was a PhD student in your group. She joined the group in 2017, ' \
    'got her MS in 2019 and PhD in 2023. She is good at math and likes good food/dancing.'

txt_dct = {
    'Joy': Joy,
    }

link_dct = {
    'Joy': 'Joy_Zhang',
    'Yalin': 'Yalin_Li',
    }


def __easter_eggs():
    cont = True
    
    while cont == True:
        name = input('Put something here... for fun: ')
        search = True
        
        if name in txt_dct.keys():
            print(f'\nShort description for {name}:')
            print(txt_dct[name])
            search = False
        if name in link_dct.keys(): 
            print(f'\nDo you want to open the webpage for {name}?')
            if input('[y]/[n]: ') in ('n', 'N', 'no', 'No', 'NO'):
                return
            webbrowser.open(prefix+link_dct[name]+suffix)
            search = False
    
        if search == True:
            webbrowser.open(f'https://www.google.com/search?d&q={name}')
        
        print('\nDo you want to continue?')
        if input('[y]/[n]: ') in ('n', 'N', 'no', 'No', 'NO'):
            return
        print('\n\n')