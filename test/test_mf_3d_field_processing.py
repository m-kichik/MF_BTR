import os
import shutil
import sys
import unittest

sys.path.extend(["../build/", "build"])
from mf_3d import MagneticField

from preparations import prepare_test_files_get_field


class TestMF3DGetField(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_dir_path = "tmp"

        if not os.path.exists(self.tmp_dir_path):
            os.mkdir(self.tmp_dir_path)

        prepare_test_files_get_field(self.tmp_dir_path)

        return super().setUp()

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_dir_path)

        return super().tearDown()
    
    def test_constant_cube(self):
        mf = MagneticField(
                os.path.join(self.tmp_dir_path, "constant_cube_field.txt")
            )
        self.assertEqual(mf.get_field(0.0, 0.0, 1.0), [42.0, 13.0, 121.0])
        self.assertEqual(mf.get_field(0.5, 0.5, 0.5), [42.0, 13.0, 121.0])
        self.assertEqual(mf.get_field(1.5, -1.5, 0.5), [42.0, 13.0, 121.0])
        self.assertEqual(mf.get_field(15, 15, 15), [42.0, 13.0, 121.0])
        self.assertEqual(mf.get_field(0.00000000001, 0.0, 0.0), [42.0, 13.0, 121.0])
        self.assertEqual(mf.get_field(0.0001, -0.0001, 0.0), [42.0, 13.0, 121.0])

    def test_opposite_vectors_cube(self):
        mf = MagneticField(
                os.path.join(self.tmp_dir_path, "opposite_vectors_cube.txt")
            )
        
        with self.subTest("Test values in nodes"):
            self.assertEqual(mf.get_field(0, 0, 0), [-1, -1, -1])
            self.assertEqual(mf.get_field(0, 0, 1), [-1, -1,  1])
            self.assertEqual(mf.get_field(0, 1, 0), [-1,  1, -1])
            self.assertEqual(mf.get_field(0, 1, 1), [-1,  1,  1])
            self.assertEqual(mf.get_field(1, 0, 0), [ 1, -1, -1])
            self.assertEqual(mf.get_field(1, 0, 1), [ 1, -1,  1])
            self.assertEqual(mf.get_field(1, 1, 0), [ 1,  1, -1])
            self.assertEqual(mf.get_field(1, 1, 1), [ 1,  1,  1])

        with self.subTest("Test cube middle"):
            self.assertEqual(mf.get_field(0.5, 0.5, 0.5), [0, 0, 0])
        
        with self.subTest("Test face middle"):
            self.assertEqual(mf.get_field(0.0, 0.5, 0.5), [-1, 0, 0])
            self.assertEqual(mf.get_field(0.5, 0.0, 0.5), [0, -1, 0])
            self.assertEqual(mf.get_field(0.5, 0.5, 0.0), [0, 0, -1])
            self.assertEqual(mf.get_field(1.0, 0.5, 0.5), [1, 0, 0])
            self.assertEqual(mf.get_field(0.5, 1.0, 0.5), [0, 1, 0])
            self.assertEqual(mf.get_field(0.5, 0.5, 1.0), [0, 0, 1])

        with self.subTest("Test edge middle"):
            self.assertEqual(mf.get_field(0.0, 0.0, 0.5), [-1, -1, 0])
            self.assertEqual(mf.get_field(0.0, 0.5, 0.0), [-1, 0, -1])
            self.assertEqual(mf.get_field(0.0, 0.5, 1.0), [-1, 0, 1])
            self.assertEqual(mf.get_field(0.0, 1.0, 0.5), [-1, 1, 0])
            self.assertEqual(mf.get_field(0.5, 0.0, 0.0), [0, -1, -1])
            self.assertEqual(mf.get_field(0.5, 0.0, 1.0), [0, -1, 1])
            self.assertEqual(mf.get_field(0.5, 1.0, 0.0), [0, 1, -1])
            self.assertEqual(mf.get_field(0.5, 1.0, 1.0), [0, 1, 1])
            self.assertEqual(mf.get_field(1.0, 0.0, 0.5), [1, -1, 0])
            self.assertEqual(mf.get_field(1.0, 0.5, 0.0), [1, 0, -1])
            self.assertEqual(mf.get_field(1.0, 0.5, 1.0), [1, 0, 1])
            self.assertEqual(mf.get_field(1.0, 1.0, 0.5), [1, 1, 0])
