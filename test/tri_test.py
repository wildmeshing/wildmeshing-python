import unittest

import os

import wildmeshing as wm


class TriangulateTest(unittest.TestCase):
    def test_run(self):
        print(wm.triangulate.__doc__)

    def test_run_no_features(self):
        root_folder = os.path.join("..", "3rdparty", "data")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        mesh_path = os.path.join(dir_path, root_folder, "rocket.obj")

        wm.triangulate(mesh_path, output="tri_test_no")

    def test_run_with_features(self):
        root_folder = os.path.join("..", "3rdparty", "data")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        mesh_path = os.path.join(dir_path, root_folder, "rocket.obj")
        feature_path = os.path.join(dir_path, root_folder, "rocket.json")

        wm.triangulate(mesh_path, feature_path, output="tri_test_with", skip_eps=True)

    # def test_svg(self):
    #     root_folder = os.path.join("..", "3rdparty", "data")
    #     dir_path = os.path.dirname(os.path.realpath(__file__))
    #     svg_path = os.path.join(dir_path, root_folder, "rocket.svg")
    #     wm.triangulate_svg(svg_path)


if __name__ == '__main__':
    unittest.main()
