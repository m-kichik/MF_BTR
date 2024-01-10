import os
import shutil
import sys
import unittest

sys.path.extend(["../build/", "build"])
from mf_3d import MagneticField

from preparations import prepare_test_files_get_field


class TestMF3DGetField(unittest.TestCase):
    def setUp(self) -> None:
        print("Call get field setup")
        self.tmp_dir_path = "tmp"

        if not os.path.exists(self.tmp_dir_path):
            os.mkdir(self.tmp_dir_path)

        prepare_test_files_get_field(self.tmp_dir_path)

        return super().setUp()

    def tearDown(self) -> None:
        print("Call get field tear down")
        shutil.rmtree(self.tmp_dir_path)

        return super().tearDown()

