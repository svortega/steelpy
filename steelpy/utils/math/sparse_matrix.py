import array
import numpy as np
import scipy.sparse as sp


class IncrementalCOOMatrix(object):

    def __init__(self, shape, dtype):

        if dtype is np.int32:
            type_flag = 'i'
        elif dtype is np.int64:
            type_flag = 'l'
        elif dtype is np.float32:
            type_flag = 'f'
        elif dtype is np.float64:
            type_flag = 'd'
        else:
            raise Exception('Dtype not supported.')

        self.dtype = dtype
        self.shape = shape

        self.rows = array.array('i')
        self.cols = array.array('i')
        self.data = array.array(type_flag)

    def append(self, i, j, v):

        m, n = self.shape

        if (i >= m or j >= n):
            raise Exception('Index out of bounds')

        self.rows.append(i)
        self.cols.append(j)
        self.data.append(v)

    def tocoo(self):

        rows = np.frombuffer(self.rows, dtype=np.int32)
        cols = np.frombuffer(self.cols, dtype=np.int32)
        data = np.frombuffer(self.data, dtype=self.dtype)

        return sp.coo_matrix((data, (rows, cols)),
                             shape=self.shape)

    def __len__(self):

        return len(self.data)
#
#
def test_incremental_coo():

    shape = 10, 10

    dense = np.random.random(shape)
    mat = IncrementalCOOMatrix(shape, np.float64)

    for i in range(shape[0]):
        for j in range(shape[1]):
            mat.append(i, j, dense[i, j])

    coo = mat.tocoo()

    assert np.all(coo.todense() == sp.coo_matrix(dense).todense())
    assert coo.row.base is mat.rows
    assert coo.col.base is mat.cols
    assert coo.data.base is mat.data
#
#
class SparseMatrixBuilder():
    def __init__(self, shape, build_size_limit):
        self.sparse_matrix = scipy.sparse.csr_matrix(shape)
        self.shape = shape
        self.build_size_limit = build_size_limit
        self.data_temp = []
        self.col_indices_temp = []
        self.row_indices_temp = []


    def add(self, data, col_indices, row_indices):
        self.data_temp.append(data)
        self.col_indices_temp.append(col_indices)
        self.row_indices_temp.append(row_indices)
        if len(self.data_temp) == self.build_size_limit:
            self.sparse_matrix += scipy.sparse.csr_matrix(
                (np.concatenate(self.data_temp),
                 (np.concatenate(self.col_indices_temp),
                  np.concatenate(self.row_indices_temp))),
                shape=self.shape
            )
            self.data_temp = []
            self.col_indices_temp = []
            self.row_indices_temp = []

    def get_matrix(self):
        self.sparse_matrix += scipy.sparse.csr_matrix(
            (np.concatenate(self.data_temp),
             (np.concatenate(self.col_indices_temp),
              np.concatenate(self.row_indices_temp))),
            shape=self.shape
        )
        self.data_temp = []
        self.col_indices_temp = []
        self.row_indices_temp = []
        return self.sparse_matrix