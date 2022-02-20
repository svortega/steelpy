# Copyright (c) 2022 steelpy
#
# Python stdlib imports
from typing import NamedTuple, Dict, List, Iterable, Union
import string

# package imports
from steelpy.process.spreadsheet.xl_pyl import ExcelInt
from steelpy.process.spreadsheet.xl_pkg import ExcelExt

#
#
class Spreadsheet:
    __slots__ = ['units', '_wb', '_wb_name',
                 '_wb_write', '_cells']
    tools:dict = {}
    def __init__(self):
        """
        """
        self._cells:dict = {}        
    #
    def __setitem__(self, key, formula):
        if isinstance(formula, str) and formula[0] == '=':
            formula = formula[1:]
        else:
            formula = (formula,)
        self._cells[key] = formula
    
    def getformula(self, key):
        c = self._cells[key]
        if isinstance(c, str):
            return '=' + c
        return c[0]
    
    def __getitem__(self, key ):
        c = self._cells[key]
        if isinstance(c, str):
            return eval(c, Spreadsheet.tools, self)
        return c[0]    
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
#
#
def col_to_num(col_str):
    """ Convert base26 column string to number. """
    expn = 0
    col_num = 0
    for char in reversed(col_str):
        col_num += (ord(char) - ord('A') + 1) * 26**expn
        expn += 1
    return col_num
#
def letter2num(letters, zbase=False):
    """A = 1, C = 3 and so on. Convert spreadsheet style column enumeration to a number.
    
    Answers:
    A = 1, Z = 26, AA = 27, AZ = 52, ZZ = 702, AMJ = 1024

    >>> letter2num('A') == 1
    True
    >>> letter2num('Z') == 26
    True
    >>> letter2num('AZ') == 52
    True
    >>> letter2num('ZZ') == 702
    True
    >>> letter2num('AMJ') == 1024
    True
    >>> letter2num('AMJ', zbase=True) == 1023
    True
    >>> letter2num('A', zbase=True) == 0
    True
    """
    letters = letters.upper()
    res = 0
    weight = len(letters) - 1
    for i, c in enumerate(letters):
        res += (ord(c) - 64) * 26**(weight - i)
    if not zbase:
        return res
    return res - 1
#
#
class ColumnName(dict):
    """
    """
    def __init__(self):
        super(ColumnName, self).__init__()
        self.alphabet = string.ascii_uppercase
        self.alphabet_size = len(self.alphabet)

    def __missing__(self, column_number):
        ret = self[column_number] = self.get_column_name(column_number)
        return ret

    def get_column_name(self, column_number):
        if column_number <= self.alphabet_size:
            return self.alphabet[column_number - 1]
        else:
            return self.alphabet[int(((column_number - 1) / self.alphabet_size)) - 1] + self.alphabet[((column_number - 1) % self.alphabet_size)]
#
