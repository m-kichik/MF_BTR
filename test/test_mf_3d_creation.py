import os
import shutil
import sys
import unittest

sys.path.extend(["../build/", "build"])
from mf_3d import MagneticField

from preparations import prepare_test_files


class TestMF3DCreation(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_dir_path = "tmp"

        if not os.path.exists(self.tmp_dir_path):
            os.mkdir(self.tmp_dir_path)

        prepare_test_files(self.tmp_dir_path)

        return super().setUp()

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_dir_path)

        return super().tearDown()

    def test_normal_cases(self):
        with self.subTest("Test file with one correct line"):
            mf = MagneticField(
                os.path.join(self.tmp_dir_path, "one_norm_line_field.txt")
            )

            self.assertEqual(mf.get_field(0.0, 0.0, 0.0), [0.0, 0.0, 0.0])

        with self.subTest("Test file with zero cube field"):
            mf = MagneticField(os.path.join(self.tmp_dir_path, "zero_cube_field.txt"))
            self.assertEqual(mf.get_field(0.0, 0.0, 1.0), [0.0, 0.0, 0.0])
            self.assertEqual(mf.get_field(0.5, 0.5, 0.5), [0.0, 0.0, 0.0])

        with self.subTest("Test file with zero cube field with empty strings"):
            mf = MagneticField(
                os.path.join(
                    self.tmp_dir_path, "zero_cube_field_with_empty_strings.txt"
                )
            )
            self.assertEqual(mf.get_field(0.0, 0.0, 1.0), [0.0, 0.0, 0.0])
            self.assertEqual(mf.get_field(0.5, 0.5, 0.5), [0.0, 0.0, 0.0])

        with self.subTest("Test file with zero cube field with different spaces"):
            mf = MagneticField(
                os.path.join(
                    self.tmp_dir_path, "zero_cube_with_different_spaces_field.txt"
                )
            )
            self.assertEqual(mf.get_field(0.0, 0.0, 1.0), [0.0, 0.0, 0.0])
            self.assertEqual(mf.get_field(0.5, 0.5, 0.5), [0.0, 0.0, 0.0])

    def test_incorrect_input(self):
        with self.subTest("Test empty file"):
            with self.assertRaises(RuntimeError):
                MagneticField(os.path.join(self.tmp_dir_path, "empty_field.txt"))

        with self.subTest("Test file with header only"):
            with self.assertRaises(RuntimeError):
                MagneticField(
                    os.path.join(self.tmp_dir_path, "empty_field_with_header.txt")
                )

        with self.subTest("Test file with one incorrect line"):
            with self.assertRaises(RuntimeError):
                MagneticField(
                    os.path.join(self.tmp_dir_path, "one_incorrect_line_field.txt")
                )

        with self.subTest("Test file with zero cube field with invalid first line"):
            with self.assertRaises(RuntimeError):
                MagneticField(
                    os.path.join(
                        self.tmp_dir_path,
                        "zero_cube_field_with_first_line_incorrect.txt",
                    )
                )

        with self.subTest("Test file with last line incorrect"):
            with self.assertRaises(RuntimeError):
                MagneticField(
                    os.path.join(
                        self.tmp_dir_path,
                        "zero_cube_field_with_last_line_incorrect.txt",
                    )
                )

        with self.subTest("Test file with middle line incorrect"):
            with self.assertRaises(RuntimeError):
                MagneticField(
                    os.path.join(
                        self.tmp_dir_path,
                        "zero_cube_field_with_middle_line_incorrect.txt",
                    )
                )


if __name__ == "__main__":
    unittest.main()