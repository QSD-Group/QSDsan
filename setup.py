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
    version='1.2.2',
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
        'biosteam>=2.34.1',
        'thermosteam>=0.30.1',
        'matplotlib>=3.3.2',
        'pandas>=1.3.2',
        'SALib>=1.4.5',
        'scikit-learn',
        'scipy>=1.7.1',
        'seaborn',
        'sympy>=1.8'
    ],
    package_data=
        {'qsdsan': [
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
    platforms=['Windows', 'Mac', 'Linux'],
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
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    keywords=[
        'quantitative sustainable design',
        'sanitation',
        'resource recovery',
        'techno-economic analysis',
        'life cycle assessment'
    ],
)