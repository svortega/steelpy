#
# Copyright (c) 2009 steelpy
#
# Python stdlib imports
from __future__ import annotations
#from typing import NamedTuple

# package imports
#from .load_case import BasicLoadSQL
#from .combination import LoadCombSQL
from steelpy.utils.sqlite.utils import create_connection, create_table
#
from steelpy.ufo.load.process.combination import (get_comb_dict,
                                                  get_comb_list)
from steelpy.ufo.load.process.load_case import (get_bload_list,
                                                get_bload_dict)
from steelpy.ufo.plot.main import PlotLoad
from steelpy.utils.dataframe.main import DBframework
#
#
#
class MasterLoad:
    """ """
    __slots__ = ['_name']

    def __init__(self, name:int|str):
        """
        """
        self._name = name


    def __str__(self) -> str:
        """ """
        output = "\n"
        output += self._basic.__str__()
        output += self._combination.__str__()
        return output

    #
    # ----------------------------
    #
    def _comb_setup(self, values: list):
        """" values : [name, type, load_id, factor, title] """
        cname = values[0]
        ltype = values[1]
        lname = values[2]
        factor = float(values[3])
        try:
            title = values[4]
        except IndexError:
            title = cname
            #
        try:
            self._combination[cname]
        except KeyError:
            self._combination[cname] = title
        #
        if re.match(r"\b(basic(_|-|\s*)?(load)?)\b", ltype, re.IGNORECASE):
            self._combination[cname]._basic[lname] = factor

        elif re.match(r"\b((load)?(_|-|\s*)?comb(ination)?(_|-|\s*)?(load)?)\b", ltype, re.IGNORECASE):
            self._combination[cname]._combination[lname] = factor

        else:
            raise IOError(f"Combination load type {ltype} not available")

    def combination(self, values: None | list | tuple | dict = None,
                    df=None):
        """
        """
        if values:
            if isinstance(values, dict):
                comb = get_comb_dict(values)
                sname = comb[2]  # lname
                if isinstance(sname, (list | tuple)):
                    db = DBframework()
                    newdf = db.DataFrame(values)
                    self._combination.df = newdf
                else:
                    self._comb_setup(values=comb)
            elif isinstance(values, (list, tuple)):
                for value in values:
                    if isinstance(value, (list, tuple)):
                        comb = get_comb_list(value)
                    elif isinstance(value, dict):
                        comb = get_comb_dict(value)
                    else:
                        raise IOError(f'load format not valid')
                    #
                    self._comb_setup(values=comb)
        #
        # dataframe input
        try:
            columns = list(df.columns)
            self._combination.df = df
        except AttributeError:
            pass
        #
        return self._combination

    #
    # -----------------------------------------------
    #
    def _basic_setup(self, values: list):
        """ [name, item, item_id, load_type, values] """
        name = values[0]
        item = values[1]
        load = values[2:]
        #
        try:
            self._basic[name]
        except IOError:
            self._basic[name] = name
        #
        if re.match(r"\b(node(s)?|point(s)?)\b", item, re.IGNORECASE):
            self._basic[name].node(load)

        elif re.match(r"\b(beam(s)?)\b", item, re.IGNORECASE):
            self._basic[name].beam(load)

        else:
            raise IOError(f"Basic load type {item[1]} not available")

    def basic(self, values: None | list | tuple | dict = None,
              df=None):
        """ """
        if values:
            if isinstance(values, dict):
                bload = get_bload_dict(values)
                itype = bload[2]  # item_type: node/beam
                if isinstance(itype, (list | tuple)):
                    db = DBframework()
                    newdf = db.DataFrame(values)
                    self._basic.df = newdf
                else:
                    self._basic_setup(values=bload)

            elif isinstance(values, (list, tuple)):
                for value in values:
                    if isinstance(value, (list, tuple)):
                        bload = get_bload_list(value)
                    elif isinstance(value, dict):
                        bload = get_bload_dict(value)
                    else:
                        raise IOError(f'load format not valid')
                    #
                    self._basic_setup(values=bload)

        #
        # dataframe input
        try:
            columns = list(df.columns)
            self._basic.df = df
        except AttributeError:
            pass
            #
        return self._basic

    #
    # ----------------------------
    #
    #
    def metocean(self, condition: None | list | tuple | dict = None,
                 df=None):
        """
        design_load : max_BS
        criterion : select design load based on local (member) or global system
        """
        #
        # if isinstance(condition, dict):
        if condition:
            # cases = []
            for key, item in condition.items():
                self._hydro[key] = item

            # 1 / 0
            # if isinstance(values[0], (list, tuple)):
            #    for value in values:
            #        self._hydro[value[0]] = value[1:]
            # else:
            #    self._hydro[values[0]] = [*values[1:], 'local']
        # elif isinstance(condition, (list, tuple)):
        #    1 / 0

        #
        # dataframe input
        try:
            columns = list(df.columns)
            header = {}
            1 / 0
        except AttributeError:
            pass
            #
        return self._hydro

    #
    # ----------------------------
    #
    # def time_history(self):
    #    """
    #    """
    #    return self.th
    #
    #
    # ----------------------------
    #
    # def mass(self):
    #    """
    #    :return:
    #    """
    #    return self._mass
    #
    #
    #
    # ----------------------------
    # Process
    # ----------------------------
    #
    #
    # ----------------------------
    # Plotting
    # ----------------------------
    #
    # @property
    def plot(self, figsize: tuple = (10, 10)):
        """ """
        return PlotLoad(cls=self, figsize=figsize)
#
#