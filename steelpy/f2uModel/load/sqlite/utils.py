

def get_load_data(conn, load_name:int|str, load_level: str):
    """ """
    cur = conn.cursor()
    if isinstance(load_name, str):
        cur.execute(f"SELECT * FROM tb_Load \
                      WHERE name = '{load_name}' AND level = '{load_level}' ;")
    else:
        cur.execute(f"SELECT * FROM tb_Load \
                      WHERE name = {load_name} AND level = '{load_level}' ;")
    loads = cur.fetchone()
    return loads