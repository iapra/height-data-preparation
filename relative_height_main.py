import os
import time

from config import LAS_INPUT_FOLDER, \
    POINTS_CSV, POLY_CSV, FILE_OUT_Z, FILE_OUT_SQL, TABLE

from relative_height import txt_SQL_cmd, get_new_z, \
    get_new_las, write_new_las


def relative_height_main(get_sql_cmd, get_intermediate_new_z, get_new_las_from_txt, \
    classification, get_new_las_directly = False):

    if get_sql_cmd == True:
        txt_SQL_cmd(LAS_INPUT_FOLDER, FILE_OUT_SQL, TABLE, classification)

    if get_intermediate_new_z == True:
        try:
            assert os.path.exists(POLY_CSV)
            assert os.path.exists(POINTS_CSV)
        except AssertionError:
            print("The csv fils containing roof polygons and points either does not exist, or not at the correct path.")
        if os.path.exists(FILE_OUT_Z):
            print(f'{FILE_OUT_Z} already exists.')
            return
            # print(f'Deleted {FILE_OUT_Z}. Writing a new one, same directory.')
            # os.remove(FILE_OUT_Z)
        get_new_z(FILE_OUT_Z, POLY_CSV, POINTS_CSV, txt_output_bool=True)

    if get_new_las_from_txt == True:
        try:
            assert os.path.exists(FILE_OUT_Z)
        except AssertionError:
            print("New z values first need to be computed. \
            If they have already been, checkout your output file and path.")
        get_new_las(LAS_INPUT_FOLDER, FILE_OUT_Z, classification)

    if get_new_las_directly == True:
        dict_points = get_new_z(FILE_OUT_Z, POLY_CSV, POINTS_CSV, txt_output_bool=False)
        write_new_las(LAS_INPUT_FOLDER, dict_points, classification)