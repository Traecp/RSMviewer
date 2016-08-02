#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Automatically import "autonomous" modules
__all__ = ['mca_spec', 'Bruker']

# # Modules depending on external (non-trivial) packages
# __others__ = []

# Set to True if you want to import all previous modules directly
importAll = True

if importAll:
    for pkg in __all__:
        __import__(__name__ + '.' + pkg)
