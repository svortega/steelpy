#
# Copyright (c) 2009 steelpy
#

# Python stdlib imports
from __future__ import annotations
import sqlite3 as sqlite3
#from sqlite3 import Error
#from datetime import datetime
#

# package imports
#


#
def create_connection(db_file):
    """ create a database connection to a SQLite database
    """
    try:
        conn = sqlite3.connect(db_file, check_same_thread=False)
        #print('sqlite version {:} {:}'.format(sqlite3.version, sqlite3.sqlite_version_info))
        return conn
    except sqlite3.Error as e:
        raise RuntimeError(e)
    #finally:
    #    conn.close()
    #return None
#
def create_table(conn, create_table_sql):
    """ create a table from the create_table_sql statement
    :param conn: Connection object
    :param create_table_sql: a CREATE TABLE statement
    :return:
    """
    #with closing(conn.cursor()) as cursor:
    #    cursor.execute(create_table_sql)
    #
    try:
        c = conn.cursor()
        c.execute(create_table_sql)
    except sqlite3.Error as e:
        raise RuntimeError(e)
#
#
#
# ------------------------------------------------------------------
#
# ------------------------------------------------------------------
# 
#