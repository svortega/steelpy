# 
# Copyright (c) 2009-2020 fem2ufo
# 


# Python stdlib imports
from itertools import chain
from collections import Counter
from collections import defaultdict
import functools
#import pickle
from typing import NamedTuple, Dict, Union, Tuple, List
import time
#

# package imports
#from fem2ufo.process.control import euclid
#import fem2ufo.process.common.operations.printing as printing
#import fem2ufo.process.common.operations.headers as headers
#from fem2ufo.f2u_model.femodel.concept.operations import find_item_fe_number_by_name

#
#
def get_nodes_connected(elements, nodes):
    """
    """
    #nodes = elements._f2u_nodes
    #print('-->')
    #
    degree = get_node_degree(elements)
    first_degree = next(iter(degree.values()))
    #
    #node_ends ={}
    connectivity = list(map(list, zip(*elements._connectivity))) 
    nodes_conn = {}
    _node_labels = []
    for key in nodes.keys():
        _node_labels.append(key)
        node_ends = [[n for n,x in enumerate(col) if x==key]
                          for col in connectivity]
        _items = []
        _memb = []
        for _col_no, _col in enumerate(node_ends):
            for _number, column in enumerate(connectivity):
                if _number == _col_no:
                    continue
                for _row in _col:
                    _items.append(column[_row])
                    _memb.append(elements._labels[_row])
        nodes_conn[key] = [_memb, _items]
    #
    levels = []
    levels.append([first_degree[0]])
    del_list_inplace(_node_labels, levels[-1])
    #_node_labels.remove(first_degree[0])
    levels.append( nodes_conn[first_degree[0]][1])
    del_list_inplace(_node_labels, levels[-1])
    rem = [levels[0][0], *levels[1]]
    # FIXME: infinite loop danger
    while _node_labels:
        rem = get_level(levels, nodes_conn, rem)
        del_list_inplace(_node_labels, levels[-1])
    #
    _node_number = 0
    _memb_number = 0
    _memb_rem = []
    for level in levels:
        # renumber node
        for _node_name in sorted(level):
            #_node_number += 1
            _index = nodes._labels.index(_node_name)
            nodes._number[_index] = _node_number
            _node_number += 1
            # renumber member
            for _member_name in sorted(nodes_conn[_node_name][0]):
                if _member_name in _memb_rem:
                    continue
                _memb_rem.append(_member_name)
                _memb_number += 1
                elements[_member_name].number = _memb_number
    #print('end')
#
#
def get_node_degree(elements):
    """
    """
    flat_nodes = list(chain.from_iterable(elements._connectivity))
    shared_nodes = Counter(flat_nodes)
    degree = defaultdict(list)
    for key, value in sorted(shared_nodes.items()):
        degree.setdefault(value, []).append(key)
    return degree
#
#
def get_level(levels, nodes_conn, rem):
    """
    """
    cases = {}
    rem2 = []
    for item in levels[-1]:
        cases[item] = set(nodes_conn[item][1]) - set(rem)
        rem2.extend(cases[item])
    rem2 = list(set(rem2))
    levels.append(rem2)
    return rem + rem2
#
def del_list_inplace(lst, id_to_del):
    for item in set(id_to_del):
        while item in lst:
            lst.remove(item)
#
#
#
@functools.lru_cache(maxsize=2048) 
def find_node_data(word_in: str) -> str:
    """
    Identify beam data from user
    """
    _key: Dict = {"number": r"\b(number|mesh)\b",
                  "name": r"\b(name|label)\b",
                  "elements": r"\b(element|member|item(s)?)\b",
                  "group": r"\b(group|set)(s)?\b",
                  "boundary": r"\b(boundar(y|ies))\b",
                  "z": r"\b(z)\b",
                  "y": r"\b(y)\b",
                  "x": r"\b(x)\b",
                  "coordinates": r"coordinates"}
    
    _match = common.find_keyword(word_in, _key)
    return _match

@functools.lru_cache(maxsize=2048) 
def find_element_data(word_in: str) -> str:
    """
    Identify beam data from user
    """
    _key: Dict = {"number": r"\b(number|mesh)\b",
                  "name": r"\b(name|label)\b",
                  "connectivity": r"\b(connectivity|node(s)?|joint(s)?)\b",
                  "material": r"\b(material)\b",
                  "section": r"\b(geometry|section|shape)\b",
                  "hinges": r"\b(hinge(s)?)\b",
                  "group": r"\b(group|set)(s)?\b",
                  # "boundary": r"\b(boundar(y|ies))\b",
                  "type": r"type"}
    
    _match = common.find_keyword(word_in, _key)
    return _match

#
def get_args(args, items, item_class, item_type):
    """
    """
    if len(args) == 1:
        for _arg in args[0]:
            if type(_arg) in [str, int, float]:
                try:
                    items[_arg]
                    raise Exception('{:} {:} already exist'
                                    .format(_arg, item_type))
                except KeyError:
                    items[_arg] = item_class(*args[0])
                break
            else:
                try:
                    items[_arg[0]]
                    raise Exception('{:} {:} already exist'
                                    .format(_arg[0], item_type))
                except KeyError:
                    items[_arg[0]] = item_class(*_arg)
    else:
        try:
            items[args[0]]
            raise Exception('{:} {:} already exist'
                            .format(args[0], item_type))
        except KeyError:
            items[args[0]] = item_class(*args)
#
#
#
#class f2uElements(NamedTuple):
#    """ Cartesian coordinate system"""
#    name:Union[str, int]
#    index:int
#    type:str
#    connectivity: List[int]
#    length:float
#    DoF:List[int]
#    beta:float
#    unit_vector: List[float]
#    material: Tuple
#    section:Tuple
#    #nodes: Dict[Union[str, int], Tuple]
##
#def dump_f2u_mesh(mesh):
#    """
#    """
#    start_time = time.time()
#    nodes = {key:item for key, item in mesh.nodes.items()}
#    material = mesh.materials.get_material()
#    sections = mesh.sections.get_properties()
#    boundaries = {key:item for key, item in mesh.boundaries.node.items()}
#    free_nodes = mesh.elements.get_free_nodes()
#    #
#    elements = {}
#    for key, memb in mesh._elements.items():
#        elements[key] = f2uElements(name=memb.name, index=memb.index, 
#                                    type= memb.type,
#                                    connectivity=memb.connectivity,
#                                    length=memb.length_node2node, 
#                                    DoF=memb.DoF, beta=memb.beta,
#                                    unit_vector=memb.unit_vector,
#                                    #K=memb.Kmatrix,
#                                    material=memb.material,
#                                    section=memb.section)
#    #
#    #  
#    file = open( "mesh.f2u", "wb" )
#    pickle.dump( nodes, file )
#    pickle.dump( boundaries, file )
#    pickle.dump( elements, file )
#    pickle.dump( material, file )
#    pickle.dump( sections, file )
#    pickle.dump( free_nodes, file )
#    file.close()
#    #
#    end_time = time.time()
#    uptime = end_time - start_time
#    print("** Writting File Process Time: {:1.4e} sec".format(uptime))
#
#
    


