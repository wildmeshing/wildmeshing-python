import unittest

import os

import wildmeshing as wm


class TriangulateTest(unittest.TestCase):
    def test_run(self):
        print(wm.triangulate.__doc__)
    
    def test_run_no_features(self):
        root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "3rdparty", "data")
        mesh_path = os.path.join(root, "rocket.obj")

        wm.triangulate(mesh_path, output="tri_test_no")

    def test_run_with_features(self):
        root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "3rdparty", "data")
        mesh_path = os.path.join(root, "rocket.obj")
        feature_path = os.path.join(root, "rocket.json")

        wm.triangulate(mesh_path, feature_path, output="tri_test_with", skip_eps=True)


if __name__ == '__main__':
    unittest.main()
