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

__all__ = ('test_auxiliary_unit_names_match_attributes',)


def test_auxiliary_unit_names_match_attributes():
    '''
    For every `SanUnit` subclass declaring a non-empty `auxiliary_unit_names`,
    each declared name must be assigned via `self.<name> = ...` somewhere in
    the class's methods (e.g., `__init__`, `_setup`, or other lifecycle methods)
    Otherwise BioSTEAM's `auxiliary_units` silently drops it (`getattr(self, name, None)`
    returns `None` with no error) without raising.
    '''
    import ast, inspect
    from textwrap import dedent
    import qsdsan as qs

    modules = [m for m in (getattr(qs, 'unit_operations', None),
                            getattr(qs, 'sanunits', None)) if m is not None]
    seen = {}
    for module in modules:
        for name in dir(module):
            obj = getattr(module, name)
            if isinstance(obj, type) and issubclass(obj, qs.SanUnit):
                seen[obj] = obj.__name__

    def assigned_attrs(cls):
        attrs = set()
        for klass in cls.__mro__:
            for method_name, member in klass.__dict__.items():
                if not callable(member):
                    continue
                try:
                    source = inspect.getsource(member)
                except (OSError, TypeError):
                    continue
                try:
                    tree = ast.parse(dedent(source))
                except (SyntaxError, IndentationError):
                    continue
                for node in ast.walk(tree):
                    # Direct assignment: self.attr = ...
                    if (isinstance(node, ast.Attribute) and isinstance(node.ctx, ast.Store)
                            and isinstance(node.value, ast.Name) and node.value.id == 'self'):
                        attrs.add(node.attr)
                    # Auxiliary unit registration: self.auxiliary('attr', ...)
                    elif (isinstance(node, ast.Call)
                            and isinstance(node.func, ast.Attribute)
                            and node.func.attr == 'auxiliary'
                            and isinstance(node.func.value, ast.Name)
                            and node.func.value.id == 'self'
                            and node.args
                            and isinstance(node.args[0], ast.Constant)):
                        attrs.add(node.args[0].value)
        return attrs

    # (class name, attr name): accepted because BioSTEAM's `auxiliary_units`
    # property does `unit = getattr(self, name, None); if unit is None: continue`
    # -- an unassigned declared name contributes zero cost/utility/construction,
    # not a silently wrong one, so this is a dead declaration, not a live bug.
    accepted_unused = {
        ('MESHDistillation', 'vacuum_system'), # never assigned anywhere in
        # BioSTEAM's MESHDistillation/its MRO; confirmed no _cost_vacuum()
        # method exists to assign it either.
        }

    bad = []
    for cls, name in seen.items():
        aux_names = getattr(cls, 'auxiliary_unit_names', ())
        if not aux_names: continue
        missing = [n for n in aux_names if n not in assigned_attrs(cls)
                   and (name, n) not in accepted_unused]
        if missing: bad.append((name, missing))

    assert not bad, (
        'these unit classes declare `auxiliary_unit_names` entries never '
        f'assigned as `self.<name> = ...` in any method: {bad}'
        )


if __name__ == '__main__':
    test_auxiliary_unit_names_match_attributes()
