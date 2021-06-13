#
# Copyright (c) 2009-2021 fem2ufo
#

# Python stdlib imports
from collections import defaultdict
from typing import NamedTuple, Dict, List, Iterable, Union
#
#
# package imports
#from steelpy.f2uModel.sql.read_model.results import 
from steelpy.f2uModel.results.sqlite.operation.process_sql import create_connection, create_table
from steelpy.f2uModel.results.operations.output import Results, print_deflections, print_member_forces
#
#
# --------------------------
#
#
class ResultSQL:
    
    __slots__ = ['_system', 'db_file', '_labels']

    def __init__(self, db_system:str,
                 db_file:Union[str,None]):
        """
        """
        self.db_file = db_file
        #
        self._create_table_node_results()
        self._create_table_element_results()
        #self._create_table_element()
    #
    #
    @property
    def node_displacement(self):
        """
        [node, load_tile, system, x, y, z, rx, ry, rz]
        """
        conn = create_connection(self.db_file)
        with conn:
            disp = self._get_node_displacements(conn)
        return disp
    #
    @node_displacement.setter
    def node_displacement(self, ndisp:List):
        """
        ndisp = [node, load_title, system, x, y, z, rx, ry, rz]
        """
        conn = create_connection(self.db_file)
        loads = self._get_load_items(conn)
        loads = {item[2]:item[0] for item in loads}
        ndisp = [[item[0], loads[item[1]], *item[2:]]
                 for item in ndisp]
        with conn:
            self.populate_node_displacements(conn, ndisp)
    #
    def print_node_displacement(self):
        """
        """
        print_deflections(self.node_displacement)
    #
    @property
    def element_force(self):
        """
        [element, load_title, node, position, system,  fx, fy, fz, mx, my, mz]
        """
        conn = create_connection(self.db_file)
        with conn:
            forces = self._get_element_forces(conn)
        return forces

    @element_force.setter
    def element_force(self, mforce:List):
        """
        mforce = [element, load_title, node, position, system,  fx, fy, fz, mx, my, mz]
        """
        conn = create_connection(self.db_file)
        loads = self._get_load_items(conn)
        loads = {item[2]:item[0] for item in loads}
        mforce = [[item[0], loads[item[1]], *item[2:]]
                  for item in mforce]        
        with conn:
            self.populate_element_force(conn, mforce)
        #print("-->")
    #
    def print_element_forces(self):
        """ """
        print_member_forces(self.element_force)
    #
    def _create_table_node_results(self) -> None:
        """ """
        _table_disp = "CREATE TABLE IF NOT EXISTS tb_ResNodeDisp(\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        node_name INTEGER NOT NULL REFERENCES tb_Nodes(name),\
                        load_number INTEGER REFERENCES tb_Load(number),\
                        system TEXT NOT NULL,\
                        x DECIMAL,\
                        y DECIMAL,\
                        z DECIMAL,\
                        rx DECIMAL,\
                        ry DECIMAL,\
                        rz DECIMAL);"
        #
        _table_reactions = "CREATE TABLE IF NOT EXISTS tb_ResNodeReactions(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            node_name INTEGER NOT NULL REFERENCES tb_Nodes(name),\
                            load_number INTEGER REFERENCES tb_Load(number),\
                            system TEXT NOT NULL,\
                            Fx DECIMAL,\
                            Fy DECIMAL,\
                            Fz DECIMAL,\
                            Mx DECIMAL,\
                            My DECIMAL,\
                            Mz DECIMAL);"
        #
        conn = create_connection(self.db_file)
        create_table(conn, _table_disp)
        create_table(conn, _table_reactions)
    #         
    def _create_table_element_results(self) -> None:
        """ """
        _table_elem_deflection = "CREATE TABLE IF NOT EXISTS tb_ResElemDeflection(\
                                    number INTEGER PRIMARY KEY NOT NULL,\
                                    element_name INTEGER NOT NULL REFERENCES tb_Elements(name),\
                                    load_number INTEGER REFERENCES tb_Load(number),\
                                    node_name INTEGER REFERENCES tb_Nodes(name),\
                                    position DECIMAL NOT NULL,\
                                    system TEXT NOT NULL,\
                                    x DECIMAL,\
                                    y DECIMAL,\
                                    z DECIMAL,\
                                    rx DECIMAL,\
                                    ry DECIMAL,\
                                    rz DECIMAL);"
        #
        _table_elm_forces = "CREATE TABLE IF NOT EXISTS tb_ResElemForce(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            element_name INTEGER NOT NULL REFERENCES tb_Elements(name),\
                            load_number INTEGER REFERENCES tb_Load(number),\
                            node_name INTEGER REFERENCES tb_Nodes(name),\
                            position DECIMAL NOT NULL,\
                            system TEXT NOT NULL,\
                            fx DECIMAL NOT NULL,\
                            fy DECIMAL NOT NULL,\
                            fz DECIMAL NOT NULL,\
                            mx DECIMAL NOT NULL,\
                            my DECIMAL NOT NULL,\
                            mz DECIMAL NOT NULL);"
        #
        _table_elm_stress = "CREATE TABLE IF NOT EXISTS tb_ResElemStress(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            element_name INTEGER NOT NULL REFERENCES tb_Elements(name),\
                            load_number INTEGER REFERENCES tb_Load(number),\
                            hot_spot INTEGER NOT NULL,\
                            sigma_x DECIMAL NOT NULL,\
                            sigma_y DECIMAL NOT NULL,\
                            sigma_z DECIMAL NOT NULL,\
                            tau_x DECIMAL NOT NULL,\
                            tau_y DECIMAL NOT NULL,\
                            tau_z DECIMAL NOT NULL);"
        #
        _table_elm_strain = "CREATE TABLE IF NOT EXISTS tb_ResElemStrain(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            hot_spot INTEGER NOT NULL);"
        #
        _table_elm_design = "CREATE TABLE IF NOT EXISTS tb_ElemDesign(\
                            number INTEGER PRIMARY KEY NOT NULL,\
                            position INTEGER NOT NULL);"
        #
        conn = create_connection(self.db_file)
        create_table(conn, _table_elem_deflection)
        create_table(conn, _table_elm_forces)
        create_table(conn, _table_elm_stress)
        create_table(conn, _table_elm_strain)
        create_table(conn, _table_elm_design)
    #
    #
    def populate_node_displacements(self, conn, value):
        """
        """
        sql = 'INSERT INTO tb_ResNodeDisp(node_name, load_number, system,\
                                          x, y, z, rx, ry, rz)\
                                          VALUES(?,?,?,?,?,?,?,?,?)'
        cur = conn.cursor()
        cur.executemany(sql, value)
    #
    def _get_node_displacements(self, conn):
        """ """
        cur = conn.cursor()
        cur.execute("SELECT name, number FROM tb_Nodes")
        n_nodes = cur.fetchall()
        n_nodes = {item[0]:item[1]-1 for item in n_nodes}
        #
        basic = self._get_displacement(conn, "basic", n_nodes)
        comb = self._get_displacement(conn, "combination", n_nodes)
        displacement = {"basic": basic, "combination": comb}
        return displacement
    #
    def _get_displacement(self, conn, load_type:str, n_nodes:Dict):
        """ """
        cur = conn.cursor()
        # Basic Load
        cur.execute("SELECT number, name, title FROM tb_Load\
                    WHERE type = '{:}'".format(load_type))
        bloads = cur.fetchall()
        #
        basic = {}
        for bl in bloads:
            cur.execute ("SELECT tb_ResNodeDisp.*\
                        FROM tb_ResNodeDisp\
                        WHERE tb_ResNodeDisp.load_number = {:}\
                        ORDER BY tb_ResNodeDisp.node_name ASC;".format(bl[0]))
            rows = cur.fetchall()
            # organize node data according to node index (number -1)
            node_res = [0] * len(n_nodes)
            for row in rows:
                index = n_nodes[row[1]]
                node_res[index] = [row[1], *row[4:]]
                #node_res.append([row[1], *row[5:]])
            basic[bl[2]] = Results(name=bl[1], number=bl[0], title=bl[2],
                                   load_type="basic", items=node_res)
        #print('-->')
        return basic
    #
    #
    def populate_element_force(self, conn, value:List):
        """ """
        sql = 'INSERT INTO tb_ResElemForce (element_name, load_number, node_name,\
                                            position, system, fx, fy, fz, mx, my, mz)\
                                            VALUES(?,?,?,?,?,?,?,?,?,?,?)'
        cur = conn.cursor()
        cur.executemany(sql, value)
        conn.commit()
        #return cur.lastrowid
    #
    def _get_element_forces(self, conn):
        """ """
        memf_basic = self._get_force(conn, "basic")
        memf_comb = self._get_force(conn, "combination")
        #print("--")
        return {"basic":memf_basic, "combination":memf_comb}
    #
    def _get_force(self, conn, load_type:str):
        """ """
        cur = conn.cursor()
        cur.execute("SELECT number, name, title FROM tb_Load \
                     WHERE type = '{:}'".format(load_type))
        loads = cur.fetchall()
        #
        elem_force = {}
        for item in loads:
            member_load = self._get_forces_items(cur, system="local", load_number=item[0])
            member_load.extend(self._get_forces_items(cur, system="global", load_number=item[0]))
            #print("--")
            elem_force[item[2]] = Results(name=item[1], number=item[0], title=item[2],
                                          load_type="basic", items = member_load)
        return elem_force
    #
    #
    def _get_forces_items(self, cur, system: str, load_number: int):
        """ """
        cur.execute ("SELECT tb_ResElemForce.element_name, tb_ResElemForce.node_name,\
                    tb_ResElemForce.position, tb_ResElemForce.system,\
                    tb_ResElemForce.fx, tb_ResElemForce.fy, tb_ResElemForce.fz,\
                    tb_ResElemForce.mx, tb_ResElemForce.my, tb_ResElemForce.mz\
                    FROM tb_ResElemForce\
                    WHERE tb_ResElemForce.system = '{:}'\
                    AND tb_ResElemForce.load_number = {:}\
                    ORDER BY tb_ResElemForce.element_name ASC;".format (system, load_number))
        rows = cur.fetchall ()
        return rows
    #
    def _get_load_items(self, conn):
        """ """
        cur = conn.cursor()
        cur.execute("SELECT * FROM tb_Load")
        loads = cur.fetchall()
        return loads
#
#
#
def get_forces_combXXX(cur, system, load_number):
    """ """
    cur.execute("SELECT tb_Elements.name, tb_ResElements.pos_x,\
                tb_ResElemForce.fx, tb_ResElemForce.fy, tb_ResElemForce.fz,\
                tb_ResElemForce.mx, tb_ResElemForce.my, tb_ResElemForce.mz\
                FROM tb_ResElements, tb_ResElemForce, tb_Elements\
                WHERE tb_ResElements.force = tb_ResElemForce.number\
                AND tb_ResElements.element_name = tb_Elements.name\
                AND tb_ResElements.system = '{:}'\
                AND tb_ResElements.lc_number = {:};"
                  .format(system, load_number))
    rows = cur.fetchall()
    return iter_items(rows)
#
#