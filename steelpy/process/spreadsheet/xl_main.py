# Copyright (c) 2021 LMSp
#
# Python stdlib imports
from dataclasses import dataclass
from pathlib import Path
import re
import os
from typing import NamedTuple, Dict, List, Iterable, Union

# package imports
from steelpy.process.spreadsheet.pylightxl.pylightxl import (readxl, readcsv, writexl, writecsv,
                                                             Database, utility_address2index)


#
#
class Spreadsheet:
    __slots__ = ['units', '_wb', '_flag']

    def __init__(self):
        """
        """
        pass
    #
    def workbook(self, wb_name:Union[str,None]=None):
        """workbook"""
        if wb_name:
            self._wb = ExcelFile(wb_name)
        else:
            self._wb = NewExcelFile()
        return self._wb
#
#
class ExcelFile:
    __slots__ = ['units', '_wb', '_ws']

    def __init__(self, wb_name: str):
        """
        """
        my_pathlib = Path(wb_name)
        self._wb = readxl(my_pathlib)
        self._ws = Sheets(self._wb)
    #
    @property
    def sheets(self):
        """ """
        return self._ws
    #
    @property
    def sheet_names(self):
        """ """
        return self._wb.ws_names
#
#
class Sheets:
    __slots__ = ['units', '_wb', '_cells']

    def __init__(self, wb: str):
        """
        """
        self._wb = wb
        #self._cells = xlCells()
    #
    def __getitem__(self, ws_name: str):
        """
        """
        s_name = self._wb.ws_names
        if ws_name in s_name:
            #return self._wb.ws(ws=ws_name)
            return xlCells(self._wb, ws_name)
        else:
            raise IOError(f"sheet {ws_name} not found")
    #
    #def __setitem__(self, ws_name: str, value) -> None:
    #    """
    #    """
    #    pass
    #
#
#
class xlCells:
    __slots__ = ['_wb', '_ws_name']

    def __init__(self, wb, ws_name: str):
        """ ws : sheet"""
        self._wb = wb
        self._ws_name = ws_name

    #
    #
    def __setitem__(self, cell_name: str, value) -> None:
        """
        """
        # if re.search(r"\:", cell_name):
        # tokens = re.split(r'[:]', cell_name)
        # row_i, col_i = utility_address2index(tokens[0])
        # row_j, col_j = utility_address2index(tokens[1])
        # for i in range(row_i, row_j+1):
        #    for j in range(col_i, col_j+1):
        #        self._wb.ws(ws=self._ws_name).update_index(row=i, col=j, val=value)
        if isinstance(value, list):
            row_i, col_j = utility_address2index(cell_name)
            columns = len(value)
            if isinstance(value[0], list):
                rows =  len(value[0])
                for j in range(columns):  # columns
                    for i in range(rows):  # rows
                        print('---', i)
                        self._wb.ws(ws=self._ws_name).update_index(row=row_i + i,
                                                                   col=col_j + j,
                                                                   val=value[j][i])
            else:
                for j in range(columns):  # columns
                    self._wb.ws(ws=self._ws_name).update_index(row=row_i,
                                                               col=col_j + j,
                                                               val=value[j])
        else:
            self._wb.ws(ws=self._ws_name).update_address(address=cell_name, val=value)

    #
    def __getitem__(self, cell_name: str):
        """
        """
        if re.search(r"\:", cell_name):
            tokens = re.split(r'[:]', cell_name)
            row_i, col_i = utility_address2index(tokens[0])
            row_j, col_j = utility_address2index(tokens[1])
            for i in range(row_i, row_j + 1):
                for j in range(col_i, col_j + 1):
                    self._wb.ws(ws=self._ws_name).update_index(row=i, col=j, val=value)
        else:
            return self._wb.ws(ws=self._ws_name).index(address=cell_name)
            # row_id, col_id = utility_address2index(cell_name)
            # return self._wb.ws(ws=self._ws_name).index(row=row_id, col=col_id)
    #
    def head(self,  row:int=1, column:Union[int,str]=1,
             title:Union[int,List]=0):
        """ """
        columns = []
        for col in self._wb.ws(ws=self._ws_name).cols:
            columns.append(col)
        #
        if isinstance(title, (list, tuple)):
            print('---')
            1/0
        else:
            headers = {col[title]: col[title+1:] for col in columns}
        #
        #print('--')
        return DataFrame(headers)
#
#
@dataclass
class DataFrame:
    __slots__ = ['_headers']

    def __init__(self, data:Union[None, Dict]=None):
        """
        """
        self._headers: Dict = {}
        if data:
            self._headers= data
    #
    def __setitem__(self, name: str, value:List) -> None:
        """
        """
        self._headers[name] = value
    #
    #
    def __getitem__(self, name: str) -> Dict :
        """
        """
        #xx = [isinstance(val, (float, int)) for val in self._headers[name]]
        #xx = any(xx)
        return self._headers[name]
    #
    # Return a string representation of self.
    def __str__(self) -> str:
        output = ""
        rows = [*([i] + j for i, j in self._headers.items())]
        gaps = [max(len(str(x)) for x in item) for item in rows]
        rows = list(zip(*rows))
        for row in rows:
            output += " ".join(f'{item:>{gaps[x]}}' for x, item in enumerate(row)) +"\n"
        return output
    #
    def drop(self):
        """ """
        pass
#
#
class DataFrameItem:
    pass
#
#
def get_data(sheet):
    """ """
    cell_obj = []
    for row in sheet.rows:
        cell_obj.append([])
        for cell in row:
            cell_obj[-1].append(cell.value)
    return cell_obj
#
#
class xlset:
    __slots__ = ['_ws']

    def __init__(self, ws):
        """ ws : sheet"""
        self._ws = ws

    #
    #
    def __setitem__(self, cell_name: str, value: float) -> None:
        """
        """
        self._ws.range(cell_name).value = value

    #
    def __getitem__(self, cell_name: int):
        """
        """
        return self._ws.range(cell_name).value
#
#
class ColumnName(dict):
    """
    """

    def __init__(self):
        import string
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
            return self.alphabet[int(((column_number - 1) / self.alphabet_size)) - 1] + self.alphabet[
                ((column_number - 1) % self.alphabet_size)]
#
#
def get_row_column(data_list: List, header: str):
    """ """
    row = [x for x, item in enumerate(data_list) if header in item]
    col = [x for x, item in enumerate(data_list[row[0]]) if item == header]
    return [row[0], col[0]]
