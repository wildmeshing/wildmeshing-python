import unittest

import wildmeshing as wm


class BendingTest(unittest.TestCase):
    def test_run(self):
        print(wm.triangulate.__doc__)


if __name__ == '__main__':
    unittest.main()
