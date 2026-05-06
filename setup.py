#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from setuptools import setup

setup(
    name='qsdsan',
    packages=['qsdsan'],
    version='1.4.4',
    license='UIUC',
    author='Quantitative Sustainable Design Group',
    author_email='quantitative.sustainable.design@gmail.com',
    description='Quantitative Sustainable Design for sanitation and resource recovery systems',
    long_description=open('README.rst', encoding='utf-8').read(),
    url='https://github.com/QSD-Group/QSDsan',
    project_urls={
        'Homepage': 'https://qsdsan.com',
        'Documentation': 'https://qsdsan.readthedocs.io',
        'Repository': 'https://github.com/QSD-Group/QSDsan',
    },
    install_requires=[
        'biosteam==2.53.10',
        'matplotlib',
        'pandas',
        'SALib',
        'scikit-learn',
        'scipy',
        'seaborn',
        'sympy',
    ],
    package_data={
        'qsdsan': [
            'data/*',
            'data/process_data/*',
            'data/sanunit_data/*',
            'data/sanunit_data/br/*',
            'data/sanunit_data/es/*',
            'data/sanunit_data/re/*',
            'equipments/*',
            'processes/*',
            'sanunits/*',
            'utils/*',
    ]},
    classifiers=[
        'License :: OSI Approved :: University of Illinois/NCSA Open Source License',
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
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.12',
    ],
    keywords=[
        'quantitative sustainable design',
        'sanitation',
        'resource recovery',
        'techno-economic analysis',
        'life cycle assessment',
    ],
)