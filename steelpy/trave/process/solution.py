# Copyright (c) 2009 steelpy

# Python stdlib imports
from __future__ import annotations

# package imports
from steelpy.trave.process.static import StaticSolver

#
# ---------------------------------
# Node displacement (u)
#
class UnSolver:
    __slots__ = ['_system', 'db_file', '_labels', '_mesh']

    def __init__(self, mesh):
        """ """
        self._mesh = mesh
        #self._load = load
        #self.db_file = db_file
        # create U node table
        #conn = create_connection(self.db_file)
        #with conn:
        #    self._create_table(conn)
        #
    #
    # ---------------------------------   
    #
    # ---------------------------------
    #
    #
    def __getitem__(self, node_name: int|str) -> tuple:
        """
        node_name : node number
        """
        try:
            self._labels.index(node_name)
            #conn = create_connection(self.db_file)
            #node = get_node(conn, node_name, self._component)
            #return node
            1 / 0
        except ValueError:
            raise IndexError('   *** node {:} does not exist'.format(node_name))   
    #
    #
    def static(self):
        """ """
        return StaticSolver(plane=self._mesh._plane)
    #
    #
    def dynamic(self):
        """ """
        1 / 0
    #
#
#
#
# ---------------------------------
#
#
#