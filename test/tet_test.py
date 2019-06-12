import unittest

import os

import wildmeshing as wm


class TetrahedralizeTest(unittest.TestCase):
    def test_doc(self):
        print(wm.tetrahedralize.__doc__)
        
    def test_run(self):
        root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "3rdparty", "data")
        mesh_path = os.path.join(root, "Octocat.obj")

        wm.tetrahedralize(mesh_path, "tet_test.msh")


if __name__ == '__main__':
    unittest.main()
