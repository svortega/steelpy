# Copyright (c) 2009 steelpy
#
# TODO: host until sure no longer needed
#
from steelpy.design.codes.main import CodeCheck
from steelpy.formulas.main import Roarks, Pilkey
from steelpy.design.main import ClampDesign as Clamp
from steelpy.metocean.main import Metocean
from steelpy.trave.main import Trave2D, Trave3D
from steelpy.trave.beam.main import Beam
from steelpy.ufo.main import UFOmodel
#
from steelpy.material.main import Material
from steelpy.sections.main import Section
from steelpy.utils.units.main import Units

from steelpy.utils.spreadsheet.main import Spreadsheet
#from steelpy.vibration.main import Vibration
#from steelpy.time_history.main import TimeHistory

# constants
__major__ = 0.  # for major interface/format changes
__minor__ = 8  # for minor interface/format changes
__release__ = 0  # for tweaks, bug-fixes, or development

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