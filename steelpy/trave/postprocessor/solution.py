# Copyright (c) 2009 steelpy

# Python stdlib imports
from __future__ import annotations

# package imports
from steelpy.utils.sqlite.utils import create_connection, create_table
from steelpy.utils.dataframe.main import DBframework


#
# ---------------------------------
# Node displacement (u)
#
class UnSQL:
    __slots__ = ['_system', 'db_file', '_labels', '_load']

    def __init__(self, load, db_file: str):
        """ """
        self._load = load
        self.db_file = db_file
        # create U node table
        conn = create_connection(self.db_file)
        with conn:
            self._create_table(conn)
    #
    # ------------------
    # SQL ops
    # ------------------
    #
    def _create_table(self, conn) -> None:
        """ """
        # Node displacement solution
        table_nodes = "CREATE TABLE IF NOT EXISTS tb_Un (\
                        number INTEGER PRIMARY KEY NOT NULL,\
                        load_name NOT NULL,\
                        load_level TEXT NOT NULL,\
                        load_system TEXT NOT NULL,\
                        node_name NOT NULL,\
                        x DECIMAL,\
                        y DECIMAL,\
                        z DECIMAL,\
                        rx DECIMAL,\
                        ry DECIMAL,\
                        rz DECIMAL);"
        #
        create_table(conn, table_nodes)
    #
    #
    # ------------------
    #
    @property
    def df(self):
        """ """
        conn = create_connection(self.db_file)
        with conn:
            Un = get_Undf(conn)       
        #
        return Un
    
    @df.setter
    def df(self, df):
        """ """
        conn = create_connection(self.db_file)
        #with conn:        
        #    cur = conn.cursor()
        #    cur.execute("SELECT tb_Nodes.name, tb_Nodes.number FROM tb_Nodes;")
        #    nodes = cur.fetchall()
        #    nodes = {item[0]:item[1] for item in nodes}
        #    #
        #    #cur.execute("SELECT tb_Elements.name, tb_Elements.number FROM tb_Elements;")
        #    #elements = cur.fetchall()
        #    #elements = {item[0]:item[1] for item in elements}            
        #    #
        #    cur.execute("SELECT tb_Load.name, tb_Load.number FROM tb_Load \
        #                 WHERE tb_Load.level = 'basic';")
        #    basic = cur.fetchall()
        #    basic = {item[0]:item[1] for item in basic}
        #
        #
        #df['load_number'] = df['name'].apply(lambda x: basic[x])
        #df['element_number'] = None # df['name'].apply(lambda x: elements[x])
        #df['node_number'] = df['node_name'].apply(lambda x: nodes[x])
        #df['system'] = 'global'
        #df['type'] = df['type'].apply(lambda x: x.lower())
        #df['type'] = 'displacement'
        #
        #df.rename(columns={"load_system": "system"}, inplace=True)        
        #
        # TODO: check this works
        # combination
        load_comb = self._load.combination()
        df_comb = load_comb.to_basic()        
        # displacments
        df = update_ndf(dfnode=df, dfcomb=df_comb)        
        #
        #
        header = ['load_name',  'load_level',
                  'load_system', #'type',
                  'node_name',
                  'x', 'y', 'z', 'rx', 'ry', 'rz']
        try: # 3d plane
            nodeconn = df[header].copy()
        except KeyError: # 2d plane
            header = ['load_name',  'load_level', 'load_system', 
                      'node_name','x', 'y', 'rz']
            nodeconn = df[header].copy()
            #nodeconn[['z', 'rx', 'ry']] = float(0)
        #
        #nodeconn.replace(to_replace=[''], value=[float(0)], inplace=True)
        #nodeconn['element_number'] = df['element_number']
        #nodeconn['title'] = df['title']
        #
        with conn:
            nodeconn.to_sql('tb_Un', conn,
                            index_label=header, 
                            if_exists='append', index=False)
        #
        #print('nodal disp results saved')
    #
#
#
# ---------------------------------
#
def get_Solution(conn):
    """ """
    #table = "SELECT tb_Load.name AS load_name, \
    #                tb_Load.title AS load_title, \
    #                tb_Load.level AS load_level, \
    #                tb_Nodes.name AS node_name, \
    #                tb_SolutionUn.* \
    #        FROM tb_Load, tb_Nodes, tb_SolutionUn \
    #        WHERE tb_SolutionUn.load_number = tb_Load.number \
    #        AND tb_SolutionUn.node_number = tb_Nodes.number "
    #
    table = "SELECT tb_Un.* FROM tb_Un"    
    #
    # Node load
    with conn:
        cur = conn.cursor()
        cur.execute(table)
        rows = cur.fetchall()
    #
    return rows
#
#
def get_Undf(conn):
    """nodes in dataframe format"""
    #
    ndata = get_Solution(conn)
    #
    #cols = ['load_name', 'load_title', 'load_level',
    #        'node_name', 'number', 
    #        'load_number', 'node_number',
    #        'load_system','load_type',
    #        'x', 'y', 'z', 'rx', 'ry', 'rz']
    #
    cols = ['number', 'load_name', 'load_level', 'load_system', 
            'node_name',
            'x', 'y', 'z', 'rx', 'ry', 'rz']    
    #
    # dataframe
    db = DBframework()
    df = db.DataFrame(data=ndata, columns=cols)
    #df['load_level'] = 'basic'
    #df = df[['load_name', 'load_level',
    #         #'load_number',  'load_title',
    #         'load_system',
    #         'node_name', #'load_type',
    #         'x', 'y', 'z', 'rx', 'ry', 'rz']]
    df.drop(labels=['number'], axis=1, inplace=True)
    return df
#
#
# ---------------------------------
#
#
#
def update_ndf(dfnode, dfcomb, 
               values:list[str] = ['x', 'y', 'z', 'rx', 'ry', 'rz']):
    """
    Update node displacements to include lcomb
    """
    db = DBframework()
    # group basic load by name
    ndgrp = dfnode.groupby('load_name')
    # get combinations with basic loads 
    #
    combgrp = dfcomb.groupby('load_name')
    for key, combfactors in combgrp:
        for row in combfactors.itertuples():
            comb = ndgrp.get_group(row.basic_load).copy()
            comb.loc[:, values] *= row.factor
            comb['load_level'] = 'combination'
            comb['load_name'] = row.load_name
            comb['load_number'] = row.load_number
            comb['load_title'] = row.load_title
            #
            try:
                dftemp = db.concat([dftemp, comb], ignore_index=True)
            except UnboundLocalError:
                dftemp = comb
    #
    try:
        #check = dftemp.groupby(['node_name', 'c']).sum().reset_index()
        dftemp = dftemp.groupby(['load_name', 'load_number','load_level',
                                 'load_title', 'load_system','node_name'],
                                  as_index=False)[values].sum()
        #test
        dfnode = db.concat([dfnode, dftemp], ignore_index=True)
    except UnboundLocalError:
        pass
    #
    return dfnode #, memb_comb
#
#