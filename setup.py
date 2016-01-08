#!/usr/bin/env ipython
# -*- coding: utf-8 -*-

"""
    Setup file.
    Source: https://docs.python.org/3.4//distutils/introduction.html

    See https://github.com/pypa/sampleproject/blob/master/setup.py for more complete example.
    http://sametmax.com/creer-un-setup-py-et-mettre-sa-bibliotheque-python-en-ligne-sur-pypi/

    To install package : python setup.py install
"""

from distutils.core import setup
setup(name='foo',
      version='1.4',
      py_modules=['sys', 'shutil', 'time', 'numpy', 'matplotlib', 'scipy'],
      )
