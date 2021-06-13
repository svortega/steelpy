#
# Copyright (c) 2009-2021 fem2ufo
# 

# Python stdlib imports
from array import array
from collections.abc import Mapping
from typing import NamedTuple, Tuple, List, Iterator, Dict, Iterable, ClassVar, Union

# package imports


#
class LoadCombinationBasic(Mapping):
    #
    def __init__(self):
        """
        """
        self._labels: array = array("I", [])
        self._title: List[str] = []
        self._number: array = array("I", [])
    #
    def __len__(self) -> int:
        return len(self._labels)
    #
    def __iter__(self):
        """
        """
        comb = set(self._labels)
        return iter(comb)
    #
    def __str__(self) -> str:
        """ """
        output = "\n"
        output += "{:}\n".format(80*"_")
        output += "\n"
        output += "{:}LOAD COMBINATIONS\n".format(25*" ")
        output += "\n"
        output += "{:}\n".format(80*".")
        output += "\n"
        for key in self._labels:
            item = self.__getitem__(key)
            # basic
            output += "Load Combination {:10s}\n".format(str(key))
            try:
                for basic_name, factor in item.basic_load.items():
                    output += "Basic       {:10s} {: 1.4e}\n".format(str(basic_name), factor)
                    #name, basic
            except TypeError:
                pass
            # comb
            try:
                for comb_name, factor in item.load_combination.items():
                    output += "Combination {:10s} {: 1.4e}\n".format(str(comb_name), factor)
            except TypeError:
                continue
            #name, comb
        #print('---')
        return output        
