import unittest
import os
import numpy as np

import wildmeshing as wm


class TriangulateTest(unittest.TestCase):
    def test_run(self):
        cubeV = np.array([[0, 0, 0],
                          [1, 0, 0],
                          [1, 1, 0],
                          [0, 1, 0],
                          [0, 0, 1],
                          [1, 0, 1],
                          [1, 1, 1],
                          [0, 1, 1]], dtype=np.float64)

        cubeF = np.array([[0, 3, 2],
                          [0, 2, 1],
                          [4, 5, 6],
                          [4, 6, 7],
                          [1, 2, 6],
                          [1, 6, 5],
                          [0, 4, 7],
                          [0, 7, 3],
                          [0, 1, 5],
                          [0, 5, 4],
                          [2, 3, 7],
                          [2, 7, 6]], dtype=np.int32)

        tetra = wm.Tetrahedralizer()
        tetra.set_mesh(cubeV, cubeF)
        tetra.tetrahedralize()
