import os


def prepare_test_files(dir_path: str):
    with open(os.path.join(dir_path, "empty_field.txt"), "w") as file:
        file.write("\n")

    with open(os.path.join(dir_path, "empty_field_with_header.txt"), "w") as file:
        file.write("MF 3D test file without field\n")

    with open(os.path.join(dir_path, "one_norm_line_field.txt"), "w") as file:
        file.writelines(
            [
                "MF 3D test file with one line field\n",
                "0 0 0 0 0 0\n",
            ]
        )

    with open(os.path.join(dir_path, "one_incorrect_line_field.txt"), "w") as file:
        file.writelines(
            [
                "MF 3D test file with one line field\n",
                "0 0 0 0 0\n",
            ]
        )

    with open(os.path.join(dir_path, "zero_cube_field.txt"), "w") as file:
        file.writelines(
            [
                "MF 3D test file with zero cube field\n",
                "0 0 0 0 0 0\n",
                "0 0 1 0 0 0\n",
                "0 1 0 0 0 0\n",
                "0 1 1 0 0 0\n",
                "1 0 0 0 0 0\n",
                "1 0 1 0 0 0\n",
                "1 1 0 0 0 0\n",
                "1 1 1 0 0 0\n",
            ]
        )

    with open(
        os.path.join(dir_path, "zero_cube_field_with_empty_strings.txt"), "w"
    ) as file:
        file.writelines(
            [
                "MF 3D test file with zero cube field\n",
                "0 0 0 0 0 0\n",
                "\n\n",
                "0 0 1 0 0 0\n",
                "0 1 0 0 0 0\n",
                "\n",
                "0 1 1 0 0 0\n",
                "1 0 0 0 0 0\n",
                "1 0 1 0 0 0\n",
                "1 1 0 0 0 0\n",
                "\n\n\n\n",
                "1 1 1 0 0 0\n",
                "\n\n\n",
            ]
        )

    with open(
        os.path.join(dir_path, "zero_cube_with_different_spaces_field.txt"), "w"
    ) as file:
        file.writelines(
            [
                "MF 3D test file with zero cube field\n",
                "0\t\t0 0 0 0 0 \n",
                "0 0 1 0 0 0\n",
                "0 1 0 0 0 0\n",
                "0 1 1\t0 0 0\n",
                "1 0 0 0 0 0\n",
                "1 0\t\t\t1 0 0 0\n",
                "1 1 0 0 0 0\t\n",
                "1 1 1 0\t\t\t0 0\n",
            ]
        )

    with open(
        os.path.join(dir_path, "zero_cube_field_with_first_line_incorrect.txt"), "w"
    ) as file:
        file.writelines(
            [
                "MF 3D test file with zero cube field\n",
                "0 0 0 0 0\n",
                "0 0 1 0 0 0\n",
                "0 1 0 0 0 0\n",
                "0 1 1 0 0 0\n",
                "1 0 0 0 0 0\n",
                "1 0 1 0 0 0\n",
                "1 1 0 0 0 0\n",
                "1 1 1 0 0 0\n",
            ]
        )

    with open(
        os.path.join(dir_path, "zero_cube_field_with_last_line_incorrect.txt"), "w"
    ) as file:
        file.writelines(
            [
                "MF 3D test file with zero cube field\n",
                "0 0 0 0 0 0\n",
                "0 0 1 0 0 0\n",
                "0 1 0 0 0 0\n",
                "0 1 1 0 0 0\n",
                "1 0 0 0 0 0\n",
                "1 0 1 0 0 0\n",
                "1 1 0 0 0 0\n",
                "1 1 1 0 0\n",
            ]
        )

    with open(
        os.path.join(dir_path, "zero_cube_field_with_middle_line_incorrect.txt"), "w"
    ) as file:
        file.writelines(
            [
                "MF 3D test file with zero cube field\n",
                "0 0 0 0 0 0\n",
                "0 0 1 0 0 0\n",
                "0 1 0 0 0 0\n",
                "0 1 1 0 0 0\n",
                "1 0 0 0 0 0\n",
                "1 0 1 0 0\n",
                "1 1 0 0 0 0\n",
                "1 1 1 0 0 0\n",
            ]
        )
