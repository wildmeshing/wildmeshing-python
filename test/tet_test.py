import unittest

import os

import wildmeshing as wm


class TetrahedralizeTest(unittest.TestCase):
    def test_doc(self):
        print(wm.tetrahedralize.__doc__)

    def test_run(self):
        root_folder = os.path.join("..", "3rdparty", "data")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        mesh_path = os.path.join(dir_path, root_folder, "small.stl")

        wm.tetrahedralize(mesh_path, "tet_test.msh", mute_log=True)
        wm.tetrahedralize(mesh_path, "tet_test.msh", mute_log=True)


if __name__ == '__main__':
    unittest.main()
