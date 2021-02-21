# Copyright (c) 2021 steelpy
#
# TODO: host until sure no longer needed
#
from steelpy.codes.main import CodeCheck
from steelpy.roarks.main import CircularRing
from steelpy.design.main import Design
#from steelpy.metocean.main import Metocean
from steelpy.trave3D.main import Trave3D
from steelpy.beam.main import Beam
from steelpy.f2uModel.main import f2uModel
from steelpy.process.units.main import Units



# constants
__major__ = 0.  # for major interface/format changes
__minor__ = 3  # for minor interface/format changes
__release__ = 3  # for tweaks, bug-fixes, or development

__version__ = '%d.%d.%d' % (__major__, __minor__, __release__)

__author__ = 'Salvador Vargas-Ortega'
__license__ = 'MIT'
__author_email__ = 'svortega@gmail.com'
#__maintainer_email__ = 'steelpy_users@googlegroups.com'
#__url__ = 'http://steelpy.readthedocs.org'
__downloadUrl__ = "https://github.com/svortega/steelpy"

# header
print('{:}'.format(52*'-'))
print('{:}   steelpy'.format(18*' '))
print('{:}   Version {:}'
      .format(14*' ',__version__))
print('{:}'.format(52*'-'))