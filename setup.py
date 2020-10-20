#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 07:07:41 2020

@author: yalinli_cabbi
"""

from setuptools import setup

setup(
    name='sanitation',
    packages=['sanitation'],
    version='0.0.5',
    license='University of Illinois/NCSA Open Source License',
    author='Sanitation Explorer Group',
    description='Module for sustainable design of non-sewered sanitation technologies',
    long_description=open('README.rst').read(),
    url="https://github.com/codesciencewater/sanitation",
    install_requires=['biosteam>=2.20.21'],
    package_data=
        {'sanitation': ['systems/*',
                        'utils/*',
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
    keywords=['WAter, Sanitation and Hygiene (WASH)', 'single-unit reinvented toilets (SURTs)', 'single-unit reinvented toilets (MURTs)', 'omni processors (OPs)', 'resource recovery', 'mass and energy balance'],
)