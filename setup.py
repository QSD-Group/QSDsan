#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 07:07:41 2020

@author: yalinli_cabbi
"""

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='sanitation',
    version='0.0.1',
    author='SanitationDevelopmentGroup',
    author_email='author@example.com',
    description='Module for sanitation unit design and simulation',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/codesciencewater/sanitation",
    packages=setuptools.find_packages(),
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
    keywords=['WAter, Sanitation and Hygiene (WASH)', 'single-unit reinvented toilets (SURT)', 'single-unit reinvented toilets (MURT)', 'resource recovery', 'mass and energy balance'],
    python_requires='>=3.6',
)