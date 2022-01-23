# Copyright (c) 2021 LMSp
#
# Python stdlib imports
from array import array
from collections.abc import Mapping
from dataclasses import dataclass
import statistics
from typing import NamedTuple, Dict, List, Iterable, Union

# package imports



#
@dataclass
class DataFrame(Mapping):
    __slots__ = ['_data']

    def __init__(self, data:Union[None, Dict, List]=None,
                 columns:Union[None, Dict, List]=None):
        """
        """
        self._data: Dict = {}
        if data:
            self._data = data
            #for key, item in data.items():
            #    self._data[key] = SeriesItem(item)
        #
        #print('--')
    #
    def __setitem__(self, name: str, value:List) -> None:
        """
        """
        self._data[name] = value
    #
    #
    def __getitem__(self, name: str) -> Dict :
        """
        """
        return SeriesItem(name, self._data[name])
        #return self._data[name]
    #
    # Return a string representation of self.
    def __str__(self) -> str:
        output = ""
        #
        #for i, j in self._data.items():
        #    new = i + " ".join(str(item) for item in j)
        #    new
        #
        rows = [*([i] + [str(item) for item in j] 
                  for i, j in self._data.items())]
        gaps = [max(len(str(x)) for x in item) for item in rows]
        rows = list(zip(*rows))
        for row in rows:
            output += " ".join(f'{item:>{gaps[x]}}' for x, item in enumerate(row)) +"\n"
        return output
    #
    def __len__(self) -> int:
        """Return the dimension of self."""
        return len(self._data)
    
    def __iter__(self):
        """
        """
        return iter(self._data)
    
    def __contains__(self, value) -> bool:
        return value in self._data    
    #
    def drop(self):
        """ """
        pass
    #
    def to_excel(self):
        """ """
        pass
    #
    def plot(self, x:List, y:List, kind:str='line'):
        """ """
        import matplotlib.pyplot as plt
        #
        y_axis = self._data[y]
        #if x :
        x_axis = self._data[x]
        #if kind.lower == 'line':
        plt.plot(x_axis, y_axis)
        #
        plt.show()         
#
#
@dataclass
class SeriesItem:
    __slots__ = ['_serie', '_name', '_dtype']

    def __init__(self, name:Union[int,float,str], data: List[float]) -> None:
        """Construct a new list object"""
        #
        self._name = name
        if all([isinstance(val, int) for val in data]):
            self._dtype = int
            self._serie: array = array('l', data[:])
        elif any([isinstance(val, str) for val in data]):
            self._dtype = str
            self._serie: List[str] = data            
        elif any([isinstance(val, float) for val in data]):
            self._dtype = float
            self._serie: array = array('d', data[:])
        else:
            self._dtype = str
            self._serie: List[str] = data
        # self._n:ClassVar = len(a) # Dimension.
    #
    def __getitem__(self, i) -> Union[float, int, str]:
        """Return the ith index of self."""
        return self._serie[i]
    #
    def __iter__(self):
        """
        """
        return iter(self._serie)
    #
    # Return a string representation of self.
    #def __str__(self) -> str:
    #    #return str(self._serie)
    #    output = [f"{self._name:}"]
    #    #if self._dtype == int:
    #    #    output.extend([f"{item: 6.0f}" for x, item in enumerate(self._serie)])
    #    #elif self._dtype == float:
    #    #    output.extend([f"{item: 1.4e}" for x, item in enumerate(self._serie)])
    #    #else:
    #    output.extend([f"{item}" for x, item in enumerate(self._serie)])
    #    gap = max([len(item) for item in output])
    #    output = [f'{item:>{gap}}' for item in output]
    #    return '\n'.join(map(str, output))
    #
    def __repr__(self):
        """Return a string representation of self."""
        args = ', '.join(repr(x) for x in self._serie)
        return '[{}]'.format(args)
    #
    def __len__(self) -> int:
        """Return the dimension of self."""
        return len(self._serie)
    
    def __add__(self, other) -> List[float]:
        """Return the sum of self and daframe object other."""
        if not isinstance(other, SeriesItem):
            raise TypeError("must be a dataframe class")
        result = [x + y for x, y in zip(self._serie, other)]
        return SeriesItem(self._name, result)   

    def __mul__(self, scalar) -> List[float]:
        """Return the product of self and numeric object alpha."""
        if not isinstance(scalar, (float, int)):
            raise TypeError("must be a scalar")
        result = [x * scalar for x in self._serie]
        return SeriesItem(self._name, result)
    
    def __truediv__(self, scalar) -> List[float]:
        """Return the product of self and numeric object alpha."""
        if not isinstance(scalar, (float, int)):
            raise TypeError("must be a scalar")
        result = [x / scalar for x in self._serie]
        return SeriesItem(self._name, result)
    #
    __rmul__ = __mul__
    #
    #
    def maxabs(self):
        """ """
        mabs = min(self._serie)
        if max(self._serie) > abs(mabs):
            mabs = max(self._serie)        
        return mabs
    #
    def max(self):
        """ """
        return max(self._serie)
    #
    def idxmax(self):
        """ """
        #max_value = max(self._serie)
        return [index for index,value in enumerate(self._serie)
                        if value== max(self._serie)]        
        #print('--')
    #
    def std(self):
        """standard deviation"""
        return statistics.stdev(self._serie)
    #
    def mean(self):
        """standard deviation"""
        return statistics.mean(self._serie)    
#