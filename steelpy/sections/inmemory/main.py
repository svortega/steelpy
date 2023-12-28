#
# Copyright (c) 2019 steelpy
#

# Python stdlib imports
from __future__ import annotations
#from collections.abc import Mapping
#

# package imports
from .angle import Angle
from .tubular import TubularIM
from .tee import Tee
from .channel import Channel
from .box import Box
from .ibeam import Ibeam
from .solid import SolidSection
#from ..process.operations import get_sect_properties
from steelpy.sections.process.utils import SectionMain
from steelpy.utils.dataframe.main import DBframework

# ---------------------------------
#
#
#
class SectionIM(SectionMain):
    __slots__ = ['_labels', '_number', '_title', '_type',
                 '_tubular', '_solid', '_ibeam', '_box',
                 '_channel', '_tee', '_angle', '_default']

    def __init__(self):
        """
        """
        #super().__init__()
        self._default: str | None = None
        self._labels: list[str|int] = []
        self._number: list[int] = []
        self._title: list[str] = []
        self._type: list = []        
        #
        self._tubular = TubularIM()
        self._solid = SolidSection()
        self._ibeam = Ibeam()
        self._box = Box()
        self._channel = Channel()
        self._tee = Tee()
        self._angle = Angle()
    #
    def __setitem__(self, shape_name: str | int,
                    properties: list[float] | dict[str, float] | str) -> None:
        """
        """
        super().__setitem__(shape_name, properties)
        shape_type = properties[0]
        self._labels.append(shape_name)
        self._type.append(shape_type)
        mnumber = next(self.get_number())
        self._number.append(mnumber)
        self._title.append(None)
    #
    #
    # def get_properties(self):
    #    """
    #    """
    #    summary = {}
    #    for key, item in self._sections.items():
    #        summary[key] = item._get_properties()
    #    return summary
    #
    #
    #@property
    #def tubular(self, shape_name:int|str):
    #    """ """
    #    return self._sections[shape_name]
    #
    #@tubular.setter
    #def tubualar(self, shape_name:int|str, properties):
    #    """ """
    #    self._sections[shape_name] = TubularBasic(name=shape_name,
    #                                              d=properties[0], t=properties[1])
    #
    # 
    @property
    def df(self):
        """ """
        db = DBframework()
        #header = ['number', 'name', 'title', 'type', 
        #          'diameter', 'wall_thickness',
        #          'height', 'web_thickness',
        #          'top_flange_width', 'top_flange_thickness',
        #          'bottom_flange_width', 'bottom_flange_thickness',
        #          'fillet_radius',
        #          'SA_inplane', 'SA_outplane',
        #          'shear_stress', 'build', 'compactness']
        #
        number = [x + 1 for x, item in enumerate(self._labels)]
        #
        datadf = {'number': number,
                  'name' : self._labels,
                  'title' : self._title,
                  'type': self._type, 
                  'diameter': [None for item in self._labels],
                  'wall_thickness': [None for item in self._labels],
                  'height': [None for item in self._labels],
                  'web_thickness': [None for item in self._labels],
                  'top_flange_width': [None for item in self._labels],
                  'top_flange_thickness': [None for item in self._labels],
                  'bottom_flange_width': [None for item in self._labels],
                  'bottom_flange_thickness': [None for item in self._labels],
                  'fillet_radius': [None for item in self._labels],
                  'SA_inplane': [1.0 for item in self._labels],
                  'SA_outplane': [1.0 for item in self._labels],
                  'shear_stress': ['maximum' for item in self._labels],
                  'build': ['welded' for item in self._labels],
                  'compactness': [None for item in self._labels]}        
        #
        #sec_df = db.DataFrame(data)
        #
        #for row in sec_df.itertuples():
        #    section = self.__getitem__(shape_name=row.name)
        #    data = section._data_df()
        #    sec_df[data.keys()].iloc[row.Index] =  data.values()
        #    #for key, item in data.items():
        #    #    sec_df[key].iloc[row.Index] = item
        #
        for idx, name in enumerate(self._labels):
            section = self.__getitem__(shape_name=name)
            data = section._data_df()
            for key, item in data.items():
                datadf[key][idx] = item
        #
        return db.DataFrame(datadf)
    
    @df.setter
    def df(self, df):
        """ """
        super().df(df)
        # Update data keeping
        self._type.extend(df['type'].tolist())
        self._labels.extend(df['name'].tolist())
        self._title.extend(df['title'].tolist())
        self._number.extend([next(self.get_number())
                             for item in self._labels])
#
