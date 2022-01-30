# Copyright (c) 2022 steelpy
#
# Python stdlib imports
from typing import NamedTuple, Dict, List, Iterable, Union

# package imports
from steelpy.process.spreadsheet.xl_pyl import ExcelInt
from steelpy.process.spreadsheet.xl_pkg import ExcelExt

#
#
class Spreadsheet:
    __slots__ = ['units', '_wb', '_wb_name',
                 '_wb_write']

    def __init__(self):
        """
        """
        pass
    #
    def read_book(self, wb_name:str):
        """read workbook"""
        # read workbook
        try:
            self._wb = ExcelInt(wb_name)
        except PermissionError:
            self._wb = ExcelExt(wb_name)
        #
        self._wb_name = wb_name
        return self._wb
    #
    def write_book(self, wb_name:Union[str,None]=None):
        """Write excel"""
        # write workbook
        if wb_name:
            try:
                self._wb_write = ExcelExt(wb_name)
            except ModuleNotFoundError:
                self._wb_write = ExcelInt(wb_name)
            else:
                1/0        
        return self._wb_write
    #
#

