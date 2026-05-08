#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import ast
from pathlib import Path


def test_sanunits_use_qsdsan_stream_facade():
    sanunits_path = Path(__file__).parents[1] / 'qsdsan' / 'sanunits'
    offenders = []

    for path in sanunits_path.glob('*.py'):
        tree = ast.parse(path.read_text(encoding='utf-8'))
        for node in ast.walk(tree):
            if not isinstance(node, ast.ImportFrom) or node.module != 'biosteam':
                continue
            if any(alias.name == 'Stream' for alias in node.names):
                offenders.append(f'{path.name}:{node.lineno}')

    assert offenders == []
