#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Sanitation Explorer: Sustainable design of non-sewered sanitation technologies
Copyright (C) 2020, Sanitation Explorer Development Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the UIUC open-source license. Please refer to
https://github.com/QSD-for-WaSH/sanitation/blob/master/LICENSE.txt
for license details.
'''

from setuptools import setup

setup(
    name='sanitation',
    packages=['sanitation'],
    version='0.0.5',
    license='University of Illinois/NCSA Open Source License',
    author='Sanitation Explorer Development Group',
    description='Module for sustainable design of non-sewered sanitation technologies',
    long_description=open('README.rst').read(),
    url="https://github.com/QSD-for-WaSH/sanitation",
    install_requires=['biosteam'],
    package_data=
        {'sanitation': [
                        'data/*',
                        'systems/*',
                        'utils/*',
                        'units/*',
                        ]},
    platforms=['Windows', 'Mac', 'Linux'],
    classifiers=['License :: OSI Approved :: University of Illinois/NCSA Open Source License',
                 'Environment :: Console',
                 'Topic :: Education',
                 'Topic :: Scientific/Engineering',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Scientific/Engineering :: Mathematics',
                 'Intended Audience :: Developers',
                 'Intended Audience :: Education',
                 'Intended Audience :: Manufacturing',
                 'Intended Audience :: Science/Research',
                 'Natural Language :: English',
                 'Operating System :: MacOS',
                 'Operating System :: Microsoft :: Windows',
                 'Operating System :: POSIX',
                 'Operating System :: POSIX :: BSD',
                 'Operating System :: POSIX :: Linux',
                 'Operating System :: Unix',
                 'Programming Language :: Python :: 3.6',
                 'Programming Language :: Python :: 3.7',
                 'Programming Language :: Python :: 3.8',
                 ],
    keywords=['WAter, Sanitation and Hygiene (WaSH)', 'single-unit reinvented toilets (SURTs)', 'single-unit reinvented toilets (MURTs)', 'omni processors (OPs)', 'resource recovery', 'mass and energy balance'],
)