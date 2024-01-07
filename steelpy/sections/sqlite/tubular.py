# 
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
#import re
#

# package imports
from steelpy.sections.utils.shape.tubular import TubularBasic
from steelpy.sections.sqlite.utils import SectionItemSQL
from steelpy.utils.sqlite.utils import create_connection
#
#from steelpy.utils.dataframe.main import DBframework
#
# ----------------------------------------
#      Standard Section Profiles
# ----------------------------------------
#
#
#
class TubularSQL(SectionItemSQL):
    """ """
    __slots__ = ['diameter', 'thickness', 'build', 'type',
                 'shear_stress', 'compactness', '_properties',
                 'FAvy', 'FAvz', 'name', 'number', 'db_file']
    
    def __init__(self, component:int, db_file:str):
        """
        Parameters
        ----------
        d : diametre
        t : wall Thickness
        Shear Stress: MAXIMUM / AVERAGE
        """
        self.db_file = db_file
        super().__init__(component=component,
                         db_file=self.db_file)
        #
    #
    #
    def __setitem__(self, shape_name: int|str, parameters: list) -> None:
        """
        parameters = []
        """
        try:
            self._labels.index(shape_name)
            raise Exception('element {:} already exist'.format(shape_name))
        except ValueError:
            #self._labels.append(shape_name)
            #
            d = parameters[0]
            t = parameters[1]
            FAvy = 1
            FAvz = 1
            shear_stress:str = 'maximum'
            build:str = 'welded'            
            compactness = None
            section = (shape_name, 
                       "Tubular",  # shape type
                       d, t,       # diameter, wall_thickess
                       None, None, # height, web_thickness
                       None, None, # top_flange_width, top_flange_thickness
                       None, None, # bottom_flange_width, bottom_flange_thickness
                       None,       # root radius
                       FAvy, FAvz,
                       shear_stress, build,
                       compactness,
                       None,)      # title
            number = self.push_section(section)
            #self._number.append(number)
    #
    def __getitem__(self, shape_name: str | int):
        """
        """
        try:
            index = self._labels.index(shape_name)
            #number = self._number[index]
        except ValueError:
            raise Exception(f" section name {shape_name} not found")
        #
        conn = create_connection(self.db_file)
        with conn:        
            row = self._pull_section(conn, shape_name)        
        #
        #
        return TubularBasic(name=row[0], 
                            diameter=row[3], thickness=row[4])    
    #
    #def geometry(self, d: Union[Units, float],
    #             t: Union[Units, float]):
    #    """
    #    Parameters
    #    ----------
    #    diameter : Diameter
    #    thickness : Wall Thickness
    #    """
    #    # Geometry
    #    self.diameter = d
    #    self.thickness = t
    #    #for key, value in kwargs.items():
    #    #    _dimensions = find_tubular_dimensions(key)
    #    #    get_dimension(self, _dimensions, value)
    #    # Geometry
    #    self.compactness = get_compactness(self)
    #
    #
    #@property
    #def d(self):
    #    """
    #    D: diameter
    #    """
    #    return self.get_item(item="diameter")
    #
    #@d.setter
    #def d(self, diameter: Union[Units, float]):
    #    """
    #    """
    #    diameter = get_sect_properties([diameter])
    #    self.update_item(item='diameter', value=diameter[0])
    #    self.push_property()
    #
    #@property
    #def t(self):
    #    """
    #    d: diametre
    #    """
    #    return self.get_item(item="wall_thickess")
    #
    #@t.setter
    #def t(self, thickness: Union[Units, float]):
    #    """
    #    """
    #    thickness = get_sect_properties([ thickness])
    #    self.update_item(item='wall_thickess', value=thickness[0])
    #    self.push_property()
    #
    #
    #def _get_section_table(self) -> tuple:
    #    """
    #    """
    #    project = (self.name, None, "Tubular",
    #               self.diameter, self.thickness,
    #               None, None,
    #               None, None,
    #               None, None,)
    #    return project
    #
    #def __getattribute__(self, attr):
    #    """
    #    Getter for myattr
    #    :param attr:
    #    :return:
    #    """
    #    if attr in self.__slots__:
    #        return self[attr]
    #    elif re.search(r"\bd\b", attr, re.IGNORECASE):
    #        return self.get_item(item="diameter")
    #    elif re.search(r"\bt(w)?\b", attr, re.IGNORECASE):
    #        return self.get_item(item="wall_thickess")
    #    else:
    #        raise AttributeError(f"Variable {attr} not found")
    #
    #def __setattr__(self, attr, value):
    #    """
    #    Setter for myattr
    #    :param attr:
    #    :return:
    #    """
    #
    #    if re.search(r"\bd\b", attr, re.IGNORECASE):
    #        diameter = get_sect_properties([value])
    #        self.update_item(item='diameter', value=diameter[0])
    #        self.push_property()
    #    elif re.search(r"t(\bw\b)?", attr, re.IGNORECASE):
    #        thickness = get_sect_properties([value])
    #        self.update_item(item='wall_thickess', value=thickness[0])
    #        self.push_property()
    #    else:
    #        raise AttributeError(f"Variable {attr} not found")
    #
    #
    #@property
    #def df(self):
    #    """ """
    #    db = DBframework()
    #    conn = create_connection(self.db_file)
    #    #with conn:
    #    df = db.read_sql_query("SELECT * FROM Section", conn)          
    #    return df 
    #
    #@df.setter
    #def df(self, df):
    #    """ """
    #    # TODO: naming seach 
    #    #mnumber = [next(self.get_number()) for _ in df.name]
    #    #
    #    df = get_sect_prop_df(df)
    #    #df.rename(columns={'d':'diameter', 'tw':'wall_thickness'},
    #    #          inplace=True)
    #    df['title'] = None 
    #    dfmat = df.reindex(columns=['name', 'title', 'type',
    #                                'diameter', 'wall_thickness'])
    #    
    #    #
    #    # push material header
    #    #dfmat = df[['name', 'type']].copy()
    #    #dfmat['title'] = None
    #    conn = create_connection(self.db_file)
    #    with conn:
    #        dfmat.to_sql('Section', conn,
    #                     index_label=['name', 'title', 'type',
    #                                  'diameter', 'wall_thickness'], 
    #                     if_exists='append', index=False)        
    #    #
    #    self._labels.extend(dfmat['name'].tolist())
    #    # TODO : fix numbering
    #    nitems = len(self._number) + 1
    #    self._number.extend([item + nitems for item in dfmat.index])
    #    #print('---')
#
#
