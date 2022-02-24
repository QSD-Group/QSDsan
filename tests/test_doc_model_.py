#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

__all__ = ('test_doc_model',)

def test_doc_model():
    from qsdsan.utils import load_example_model
    model = load_example_model(evaluate=True, file=None)


if __name__ == '__main__':
    test_doc_model()