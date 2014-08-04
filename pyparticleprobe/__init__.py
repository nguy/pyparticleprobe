"""
muPhysPy - Python Microphysical Analysis Package
==================================
Top-level package (:mod:`muphyspy`)
==================================

.. currentmodule:: muphyspy

.. autosummary::
    :toctree: generated/

"""
from .ReadFile import read_file
#from . import radar

__all__ = [s for s in dir() if not s.startswith('_')]
