import os


def prepare_test_files_creation(dir_path: str):
    prepare_correct_cases(dir_path)
    prepare_incorrect_cases(dir_path)


def prepare_test_files_get_field(dir_path: str):
    with open(os.path.join(dir_path, "constant_cube_field.txt"), "w") as file:
        file.writelines(
            [
                "MF 3D test file with zero cube field\n",
                "0 0 0 42 13 121\n",
                "0 0 1 42 13 121\n",
                "0 1 0 42 13 121\n",
                "0 1 1 42 13 121\n",
                "1 0 0 42 13 121\n",
                "1 0 1 42 13 121\n",
                "1 1 0 42 13 121\n",
                "1 1 1 42 13 121\n",
            ]
        )

    with open(os.path.join(dir_path, "opposite_vectors_cube.txt"), "w") as file:
        file.writelines(
            [
                "MF 3D test file with zero cube field\n",
                "0 0 0 -1 -1 -1\n",
                "0 0 1 -1 -1  1\n",
                "0 1 0 -1  1 -1\n",
                "0 1 1 -1  1  1\n",
                "1 0 0  1 -1 -1\n",
                "1 0 1  1 -1  1\n",
                "1 1 0  1  1 -1\n",
                "1 1 1  1  1  1\n",
            ]
        )


def prepare_correct_cases(dir_path: str):
    with open(os.path.join(dir_path, "one_norm_line_field.txt"), "w") as file:
        file.writelines(
            [
                "MF 3D test file with one line field\n",
                "0 0 0 1 1 1\n",
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
        os.path.join(dir_path, "zero_cube_field_without_hatter.txt"), "w"
    ) as file:
        file.writelines(
            [
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

    with open(os.path.join(dir_path, "zero_cube_field_with_hatter.txt"), "w") as file:
        file.writelines(
            [
                "MF 3D test file with zero cube field\n",
                "With hatter\n",
                "With a number of lines\n",
                "Just for testing\n",
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


def prepare_incorrect_cases(dir_path: str):
    with open(os.path.join(dir_path, "empty_field.txt"), "w") as file:
        file.write("\n")

    with open(os.path.join(dir_path, "empty_field_with_header.txt"), "w") as file:
        file.write("MF 3D test file without field\n")

    with open(os.path.join(dir_path, "one_incorrect_line_field.txt"), "w") as file:
        file.writelines(
            [
                "MF 3D test file with one line field\n",
                "0 0 0 0 0\n",
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

    with open(
        os.path.join(dir_path, "zero_cube_field_with_missed_line_first.txt"), "w"
    ) as file:
        file.writelines(
            [
                "MF 3D test file with zero cube field\n",
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
        os.path.join(dir_path, "zero_cube_field_with_missed_line_middle.txt"), "w"
    ) as file:
        file.writelines(
            [
                "MF 3D test file with zero cube field\n",
                "0 0 0 0 0 0\n",
                "0 0 1 0 0 0\n",
                "0 1 0 0 0 0\n",
                "0 1 1 0 0 0\n",
                "1 0 0 0 0 0\n",
                "1 1 0 0 0 0\n",
                "1 1 1 0 0 0\n",
            ]
        )

    with open(
        os.path.join(dir_path, "zero_cube_field_with_missed_line_last.txt"), "w"
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
            ]
        )
