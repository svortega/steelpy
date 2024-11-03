#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
from collections import defaultdict
import re

# package imports
import steelpy.utils.io_module.text as common
from steelpy.sections.utils.shape.ibeam import get_Isection, get_Isect_dict
from steelpy.sections.utils.shape.tubular import get_tub_section, get_TubSect_dict
from steelpy.sections.utils.shape.box import get_box_section, get_BoxSect_dict
from steelpy.sections.utils.shape.channel import get_Csection, get_Csect_dict
from steelpy.sections.utils.shape.tee import get_Tsection, get_Tsect_dict
from steelpy.sections.utils.shape.angle import get_Lsection, get_Lsect_dict
from steelpy.sections.utils.shape.solid import get_solid_section, get_SolidSect_dict

#from steelpy.sections.utils.operations import get_sect_properties
from steelpy.sections.utils.shape.tubular import TubularBasic
from steelpy.sections.utils.shape.solid import RectangleSolid, CircleSolid, TrapezoidSolid
from steelpy.sections.utils.shape.ibeam import IbeamBasic, get_Isection
from steelpy.sections.utils.shape.box import BoxBasic
from steelpy.sections.utils.shape.channel import ChannelBasic
from steelpy.sections.utils.shape.tee import TeeBasic
from steelpy.sections.utils.shape.angle import AngleBasic
#
#
from steelpy.utils.dataframe.main import DBframework

#
# ---------------------------
#
def get_section(properties:tuple|list|dict)->list:
    """ """
    section_type = None
    if isinstance(properties, dict):
        matkeys = list(properties.keys())
        for key in matkeys:
            if re.match(r"\b((sect(ion)?(_|\s*)?)?type|shape)\b", key, re.IGNORECASE):
                section_type = find_sect_type(properties[key])
                properties.pop(key)
                break
    elif isinstance(properties, (list, tuple)):
        section_type = properties.pop(0)
        section_type = find_sect_type(section_type)
    else:
        raise IOError(f' section input not recognised')
    #
    if section_type:
        match section_type:
            case 'I section':
                section = get_Isection(properties)
            case 'Tee':
                section = get_Tsection(properties)
            case 'Tubular':
                section = get_tub_section(properties)
            case 'Box':
                section = get_box_section(properties)
            case 'Channel':
                section = get_Csection(properties)
            case 'Angle':
                section = get_Lsection(properties)
            case 'Rectangle bar':
                section = get_solid_section(section_type, properties)
            case 'Trapezoid bar':
                section = get_solid_section(section_type, properties)
            case 'Circular bar':
                section = get_solid_section(section_type, properties)
            case 'General section':
                raise NotImplementedError()
    else:
        raise IOError(f'Section type: {section_type} not valid')
    #
    return section
#
# ---------------------------------
#
#
def find_sect_type(word_in:str):
    """
    """
    key = {"I section": r"\b((i((_|-|\s*)?beam|section)?|w|m|s|hp|ub|uc|he|ipe|pg|asb)(_|-|\s*)?(asymmetrical)?)\b",
           "Tee": r"\b(t(ee)?((_|-|\s*)?section)?)\b",
           "Tubular": r"\b(tub(ular)?|pipe|chs((_|-|\s*)?section)?)\b",
           "Box": r"\b(b(ox)?|rhs|shs((_|-|\s*)?section)?)\b",
           "Channel": r"\b(c(hannel)?((_|-|\s*)?section)?)\b",
           "Angle": r"\b(l|angle((_|-|\s*)?section)?)\b",
           "Rectangle bar": r"\b((solid|bar(_|-|\s*)?)?(square|rectangle)(_|-|\s*)?(bar)?)\b",
           "Trapezoid bar": r"\b((solid|bar(_|-|\s*)?)?trapezoid(al)?(_|-|\s*)?(bar)?)\b",
           "Circular bar": r"\b((solid|bar(_|-|\s*)?)?(circular|round)(_|-|\s*)?(bar)?)\b",
           "General section": r"\b(general((_|-|\s*)?section)?)\b",}
    match = common.find_keyword(word_in, key)
    return match
#
def find_parameters(word_in:str) -> str:
    """ """
    key = {"h": r"\b(h(eight)?)\b",
           "tw": r"\b(t(hickness|hk)?(_|-|\s*)?w(eb|all)?)\b",
           "bf": r"\b((b(ase)?|w(idth)?)(_|-|\s*)?f(lange)?((_|-|\s*)?t(op)?)?)\b",
           "tf": r"\b(t(hickness|hk)?(_|-|\s*)?f(lange)?((_|-|\s*)?t(op)?)?)\b",
           "bfb": r"\b((b(ase)?|w(idth)?)(_|-|\s*)?f(lange)?(_|-|\s*)?b(ottom)?)\b",
           "tfb": r"\b(t(hickness|hk)?(_|-|\s*)?f(lange)?(_|-|\s*)?b(ottom)?)\b",
           "r": r"\b(r(oot)?(_|-|\s*)?(r(adius|atio)?)?)\b",
           # tubular/circular
           "d": r"\b(d(iamet(ro|er|re)|epth)?)\b",
           # solid
           "a": r"\b(a)\b",
           "b": r"\b(b(ase)?)\b" }
    
    match = common.find_keyword(word_in, key)
    return match
#
#
def find_sect_item(word_in:str) -> str:
    """ """
    key = {"name": r"\b(id|name|load(s)?)\b",
           "type": r"\b(((sect(ion)?|shape)(_|-|\s*)?)?type)\b",
           "title": r"\b(title|comment)\b",
           "SA_inplane": r"\b(s(hear)?(_|-|\s*)?a(rea)?(_|-|\s*)?in(_|-|\s*)?plane)\b",
           "SA_outplane": r"\b(s(hear)?(_|-|\s*)?a(rea)?(_|-|\s*)?out(_|-|\s*)?plane)\b",
           "shear_stress": r"\b(shear(_|-|\s*)?stress)\b",
           "build": r"\b(build)\b",
           "compactness": r"\b(compactness)\b",}
    try:
        match = common.find_keyword(word_in, key)
        return match
    except IOError:
        return find_shape_item(word_in)
#
def find_shape_item(word_in: str) -> str:
    """ """
    # Basic shapes
    #try:
    #    return find_sect_type(word_in)
    #except IOError:
    #    pass
    try:
        return find_parameters(word_in)
    except IOError:
        pass

    #
    raise IOError('load definition')
#
#
#
def get_sect_properties(properties:list[float], steps:int=10):
    """ """
    #sect_prop = [None for _ in range(steps)]
    sect_prop = []
    for x, item in enumerate(properties):
        try:
            #sect_prop[x] = item.value
            sect_prop.append(item.value)
        except AttributeError:
            #sect_prop[x] = item
            sect_prop.append(item)
            raise IOError('units required')
    return sect_prop
#
#
def get_sect_prop_df(df):
    """ """
    for cname in df:
        try:
            df[cname] = df[cname].apply(lambda x: x.convert('metre').value)
        except AttributeError:
            pass
    #print('---')
    return df
#
#
def get_sect_df(df):
    """ """
    #
    columns = list(df.columns)
    header = {key: find_sect_item(key)
              for key in columns}
    df.rename(columns=header, inplace=True)
    df['type'] = df['type'].apply(lambda x: find_sect_type(x))
    #
    columns = list(df.columns)
    if 'SA_inplane' not in columns:
        df['SA_inplane'] = float(1.0)
    
    if 'SA_outplane' not in columns:
        df['SA_outplane'] = float(1.0)
    
    if 'shear_stress' not in columns:
        df['shear_stress'] = 'maximum'
    
    if 'build' not in columns:
        df['build'] = 'welded'
    
    if 'compactness' not in columns:
        df['compactness'] = None
    
    if 'title' not in columns:
        df['title'] = None
    #
    group = df.groupby("type")
    newgrp = [] # defaultdict(list)
    newprop = []
    for shape_type, items in group:
        section = items.to_dict(orient='split', index=False)
        match shape_type:
            case 'I section':
                # [d, tw, bf, tf, bfb, tfb, r,
                # FAvy, FAvz, shear_stress,
                # build, compactness, title]
                for load in section['data']:
                    data = dict(zip(section['columns'], load))                
                    line, prop = get_Isect_dict(data)
                    empty = [None] * 2
                    newgrp.append([data['name'], data['type'], 
                                   *empty, # d' tw
                                   *line])
                    newprop.append([data['name'], *prop[:]])                    
            case 'Tee':
                #[shape, d/h, tw, b, tw, r,
                # FAvy, FAvz, shear_stress,
                # build, compactness, title]
                for load in section['data']:
                    data = dict(zip(section['columns'], load))                
                    line, prop = get_Tsect_dict(data)
                    empty = [None] * 2
                    newgrp.append([data['name'], data['type'], 
                                   *empty, # d' tw
                                   *line[:4], *empty,
                                   *line[4:]])
                    newprop.append([data['name'], *prop[:]])
            case 'Tubular':
                # [diameter, thickness, 
                # FAvy, FAvz, shear_stress, 
                # build, compactness, title]
                for load in section['data']:
                    data = dict(zip(section['columns'], load))
                    line, prop = get_TubSect_dict(data)
                    empty = [None] * 7
                    newgrp.append([data['name'], data['type'], 
                                   *line[:2], # d' tw
                                   *empty, *line[2:]])
                    newprop.append([data['name'], *prop[:]])               
            case 'Box':
                #[shape, d/h, tw, b, tw, r,
                # FAvy, FAvz, shear_stress,
                # build, compactness, title]
                for load in section['data']:
                    data = dict(zip(section['columns'], load))                
                    line, prop = get_BoxSect_dict(data)
                    empty = [None] * 2
                    newgrp.append([data['name'], data['type'], 
                                   *empty, # d' tw
                                   *line[:4], *empty,
                                   *line[4:]])
                    newprop.append([data['name'], *prop[:]])
            case 'Channel':
                #[shape, d/h, tw, b, tw, r,
                # FAvy, FAvz, shear_stress,
                # build, compactness, title]
                for load in section['data']:
                    data = dict(zip(section['columns'], load))                
                    line, prop = get_Csect_dict(data)
                    empty = [None] * 2
                    newgrp.append([data['name'], data['type'], 
                                   *empty, # d' tw
                                   *line[:4], *empty,
                                   *line[4:]])
                    newprop.append([data['name'], *prop[:]])
            case 'Angle':
                #[d/h, tw, b, tw, r,
                # FAvy, FAvz, shear_stress,
                # build, compactness, title]
                for load in section['data']:
                    data = dict(zip(section['columns'], load))                
                    line, prop = get_Lsect_dict(data)
                    empty = [None] * 2
                    newgrp.append([data['name'], data['type'], 
                                   *empty, # d' tw
                                   *line[:4], *empty,
                                   *line[4:]])
                    newprop.append([data['name'], *prop[:]])
            case 'Rectangle bar':
                #[diameter/height, base, a, c,
                # FAvy, FAvz, shear_stress,
                # build, compactness, title]
                for load in section['data']:
                    data = dict(zip(section['columns'], load))
                    line, prop = get_SolidSect_dict(shape='Rectangle bar',
                                                    parameters=data)
                    empty = [None] * 2
                    newgrp.append([data['name'], data['type'],
                                   *empty, # d' tw
                                   line[0], None,
                                   line[1], None,
                                   line[2], None, None,
                                   *line[4:]])
                    newprop.append([data['name'], *prop[:]])
            case 'Trapezoid bar':
                #[diameter/height, base, a, c,
                # FAvy, FAvz, shear_stress,
                # build, compactness, title]
                for load in section['data']:
                    data = dict(zip(section['columns'], load))
                    line, prop = get_SolidSect_dict(shape='Trapezoid bar',
                                                    parameters=data)
                    empty = [None] * 2
                    newgrp.append([data['name'], data['type'],
                                   *empty, # d' tw
                                   line[0], None,
                                   line[1], None,
                                   line[2], None, None,
                                   *line[4:]])
                    newprop.append([data['name'], *prop[:]])
            case 'Circular bar':
                #[diameter/height, base, a, c,
                # FAvy, FAvz, shear_stress,
                # build, compactness, title]
                for load in section['data']:
                    data = dict(zip(section['columns'], load))
                    line, prop = get_SolidSect_dict(shape='Circular bar',
                                                    parameters=data)
                    empty = [None] * 2
                    newgrp.append([data['name'], data['type'],
                                   line[0], None, # d' tw
                                   *empty, *empty, *empty, None,
                                   *line[4:]])
                    newprop.append([data['name'], *prop[:]])
            case 'General section':
                raise NotImplementedError()
    #
    # Sections
    header = ['name', 'type', 
              'diameter', 'wall_thickness',
              'height', 'web_thickness',
              'top_flange_width', 'top_flange_thickness',
              'bottom_flange_width', 'bottom_flange_thickness',
              'fillet_radius', #'web_orientation', 
              'SA_inplane', 'SA_outplane',
              'shear_stress', 'build', 'compactness', 'title']
    db = DBframework()
    sectdf = db.DataFrame(data=newgrp, columns=header)
    #
    # Properties
    header = ['name', 'area', 'Zc', 'Yc',
              'Iy', 'Zey', 'Zpy', 'ry', 
              'Iz', 'Zez', 'Zpz', 'rz',
              'J', 'Cw',
              'alpha_sy', 'alpha_sz', 'density']
    propdf = db.DataFrame(data=newprop, columns=header)
    propdf.drop(columns=['alpha_sy', 'alpha_sz', 'density'])
    #
    return sectdf, propdf
#
#
# ---------------------------------
#
#
def ShapeGeometry(section_type: str | int, geometry: list):
    """ """
    section_type = find_sect_type(section_type)
    match section_type:
        case 'Tubular':
            return TubularBasic(name=geometry[0],
                                diameter=geometry[3],
                                thickness=geometry[4])
        case 'I section':
            return IbeamBasic(name=geometry[0],
                              d=geometry[5], tw=geometry[6],
                              bft=geometry[7], tft=geometry[8],
                              bfb=geometry[9], tfb=geometry[10])
        case 'Box':
            return BoxBasic(name=geometry[0],
                            d=geometry[5], tw=geometry[6],
                            b=geometry[7], tb=geometry[8])
        case 'Channel':
            return ChannelBasic(name=geometry[0],
                                d=geometry[5], tw=geometry[6],
                                b=geometry[7], tb=geometry[8])
        case 'Tee':
            return TeeBasic(name=geometry[0],
                            d=geometry[5], tw=geometry[6],
                            b=geometry[7], tb=geometry[8])
        case 'Angle':
            return AngleBasic(name=geometry[0],
                              d=geometry[5], tw=geometry[6],
                              b=geometry[7], r=0)
        case 'Circular bar':
            d = geometry[5]
            return CircleSolid(name=geometry[0], d=d)
        case 'Rectangle bar':
            d = geometry[5]
            wb = geometry[7]
            return RectangleSolid(name=geometry[0], depth=d, width=wb)
        case 'Trapezoid bar':
            d = geometry[5]
            wb = geometry[7]
            wt = geometry[9]
            c = abs(wt - wb) / 2.0
            return TrapezoidSolid(name=geometry[0], depth=d, width=wb, a=wt, c=c)
        case 'General section':
            raise NotImplementedError()
        case _:
            raise IOError(f' Section type {section_type} not recognised')
#