#!/usr/bin/env python

"""
setup.py file for aum
"""

from distutils.core import setup, Extension

cosmology_module = Extension('_cosmology',
                           sources=['src/cosmology.i', 'src/cosmology.cpp','src/gauleg.cpp','src/haloes.cpp','src/powerspectrum.cpp'],
                           swig_opts=["-c++"],
                           libraries=['gsl','gslcblas','m'],
                           )

setup (name        = 'aum',
       version     = '1.01',
       author      = "Surhud More",
       url         = "http://member.ipmu.jp/surhud.more/research",
       author_email= "surhud.more@ipmu.jp",
       description = """Cosmology module""",
       ext_modules = [cosmology_module],
       license     = ['GPL'],
       py_modules  = ["cosmology"],
       package_dir = { '':'src'},
       )

