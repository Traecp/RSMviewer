#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Tra NGUYEN THANH <thanh-tra.nguyen@esrf.fr>'
__version__ = '1.17'
__adv__ = 'setup.py'

from distutils.core import setup

# Package name
name = 'RSM_Viewer'

# Packages (subdirectories in lib/)
packages = ['']

# Modules (files in lib/)
modules = ['mca_spec', 'Bruker']

# Scripts (in scripts/)
scripts = ['RSMviewer.py']

cmdclass = {}
command_options = {}


setup(name=name,
      description='GUI tool for RSM analysis from the linear Vantec detector or from a labsource instrument (Bruker, Panalytical)',
      author='Tra NGUYEN THANH',
      author_email='thanh-tra.nguyen@esrf.fr',
      maintainer="Tra",
      maintainer_email='thanh-tra.nguyen@esrf.fr',
      url='http://esrf.fr/',
      license="CeCILL-C FREE SOFTWARE LICENSE",
      package_dir={name: 'lib'},
      packages=packages,
      py_modules=['.'.join((name, mod)) for mod in modules], # Individual modules
      scripts=['scripts/'+script for script in scripts],
      # Data - setuptools specific
      cmdclass=cmdclass,
      command_options=command_options
      )
