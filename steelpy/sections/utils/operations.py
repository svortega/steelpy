#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
import re

# package imports
import steelpy.utils.io_module.text as common
from steelpy.sections.utils.shape.ibeam import get_Isection
from steelpy.sections.utils.shape.tubular import get_tub_section
from steelpy.sections.utils.shape.box import get_box_section
from steelpy.sections.utils.shape.channel import get_Csection
from steelpy.sections.utils.shape.tee import get_Tsection
from steelpy.sections.utils.shape.angle import get_Lsection
from steelpy.sections.utils.shape.solid import get_solid_section



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
            case 'Tee section':
                section = get_Tsection(properties)
            case 'Tubular section':
                section = get_tub_section(properties)
            case 'Box section':
                section = get_box_section(properties)
            case 'Channel section':
                section = get_Csection(properties)
            case 'Angle section':
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
#
def find_sect_type(word_in:str):
    """
    """
    key = {"I section": r"\b(i((_|-|\s*)?beam|section)?|w|m|s|hp|ub|uc|he|ipe|pg|asb)\b",
           "Tee section": r"\b(t(ee)?((_|-|\s*)?section)?)\b",
           "Tubular section": r"\b(tub(ular)?|pipe|chs((_|-|\s*)?section)?)\b",
           "Box section": r"\b(b(ox)?|rhs|shs((_|-|\s*)?section)?)\b",
           "Channel section": r"\b(c(hannel)?((_|-|\s*)?section)?)\b",
           "Angle section": r"\b(l|angle((_|-|\s*)?section)?)\b",
           "Rectangle bar": r"\b((solid|bar(_|-|\s*)?)?square|rectangle(bar(_|-|\s*)?)?)\b",
           "Trapezoid bar": r"\b((solid|bar(_|-|\s*)?)?trapezoid(bar(_|-|\s*)?)?)\b",
           "Circular bar": r"\b((solid|bar(_|-|\s*)?)?circular|round(bar(_|-|\s*)?)?)\b",
           "General section": r"\b(general((_|-|\s*)?section)?)\b",}
    match = common.find_keyword(word_in, key)
    return match
#
#
# ---------------------------------
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