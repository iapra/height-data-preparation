import shutil
import geopandas as gpd
import os
import numpy as np
import scipy.spatial
import random
import numpy as np
import math

from config import DIR_IMAGES_GEOTIFF, NUMBER_OF_SPLITS, DATA_SPLIT_FOLDER, \
    PERCENTAGE_TRAIN, PERCENTAGE_VAL, PERCENTAGE_TEST

from utils import get_img_bbox, datasets_to_geojson

### ENTER WHAT YOU NEED
create_split = False
visualise_split = True
DATA_SPLIT_FOLDER_TO_VISUALIZE = DATA_SPLIT_FOLDER + "_3"
###########################

def train_val_test_split(DIR_IMAGES_GEOTIFF, NUMBER_OF_SPLITS, PERCENTAGE_TRAIN, PERCENTAGE_VAL, PERCENTAGE_TEST):
    '''
    Function that splits data input into train-val-test sets, where 3 datasets MIGHT OVERLAP.
    
    input: DIR_IMAGES_GEOTIFF = folder directory of geotif images
    NUMBER_OF_SPLITS = number of split wished
    PERCENTAGE_TRAIN = percentage of data to set as training data
    PERCENTAGE_VAL = percentage of data to set as validation data
    PERCENTAGE_TEST = percentage of data to set as test data
    
    output: a folder per split, containing 3 txt files (train, validation, test)
    '''
    # WE GATHER AND KEEP TRACK OF IMAGES' CENTERS
    dict_img = {}
    dict_img_inverted = {}
    center_x, center_y = [], []
    count = 0
    for image in os.listdir(DIR_IMAGES_GEOTIFF):
        if image[-4:] == ".tif":
            img_path = f"{DIR_IMAGES_GEOTIFF}/{image}"
            _, _, _, [ulx_deg, uly_deg, lrx_deg, lry_deg] = get_img_bbox(img_path)
            x = np.float128(ulx_deg + (lrx_deg-ulx_deg))
            y = np.float128(lry_deg + (lry_deg-uly_deg))
            # Assert that there are no duplicates (name or content)
            # Duplicates will be ignored: Only a unique image (last image encountered) is kept.
            try:
                assert image[:-4] not in dict_img
                try:
                    # TWO ASSERTIONS ARE PASSED, THE IMAGE IS KEPT
                    assert (x,y) not in dict_img_inverted
                    center_x.append(x)
                    center_y.append(y)
                    dict_img[image[:-4]] = (x, y)
                    dict_img_inverted[(x,y)] = image[:-4]
                    count += 1
                except:
                    continue
                    print(f"All images are not unique, check out your input. E.g. {dict_img_inverted[(x,y)]}.tif and {image[:-4]}.tif are similar.")
            except:
                continue
                print(f"All images do not have a unique name, check out your input. E.g. {image[:-4]}.tif exists at least twice.")
    try:
        assert count == len(dict_img)
    except:
        print("Please, check out your input images.")

    # KD-TREE FROM IMG CENTERS
    centers = np.array(list(zip(center_x,center_y)))
    X_tree = scipy.spatial.KDTree(centers)

    # WE CHOOSE RANDOM IMGS AS SPLIT-CENTERS
    split_centers = random.sample(list(dict_img), NUMBER_OF_SPLITS)
    print(f"The image split centers are {split_centers}.tif")

    # HOW MANY TRAIN/VAL/TEST SAMPLES
    train_nb = round(len(dict_img) * PERCENTAGE_TRAIN)
    val_nb = round((len(dict_img) - train_nb) * PERCENTAGE_VAL/(PERCENTAGE_VAL+PERCENTAGE_TEST))
    test_nb = (len(dict_img) - (val_nb + train_nb))
    print(f"There are {len(dict_img)} images. The split has {train_nb} training images, {val_nb} validation images, {test_nb} test images.")
    assert (train_nb + val_nb + test_nb == len(dict_img))

    # WE RETRIEVE THE 60% CLOSEST POINT 
    for i in range(NUMBER_OF_SPLITS):
        # create folder structure, if not existing
        FOLDER = DATA_SPLIT_FOLDER + '_' + str(i+1)
        if os.path.isdir(FOLDER):
            shutil.rmtree(FOLDER)
        if not os.path.isdir(FOLDER):
            os.mkdir(FOLDER)

        split_center = split_centers[i]
        split_center_x, split_center_y = dict_img[split_center]
        center = (split_center_x, split_center_y)
        distance, indices = X_tree.query(center, train_nb)

        # IN EACH SPLIT-FOLDER WE CREATE 3 TXT FILES
        open(f'{FOLDER}/train.txt', 'w').close
        open(f'{FOLDER}/val.txt', 'w').close
        open(f'{FOLDER}/test.txt', 'w').close

        # WE APPEND THE TXT FILES WITH CORRESPONDING IMAGES
        assert len(indices) == train_nb 

        ### TRAINING SET ###
        for index in indices:
            center_coords = tuple(centers[index])
            image_name = dict_img_inverted[center_coords]

            with open(f'{FOLDER}/train.txt', 'a') as f1:
                f1.write(f"{int(image_name)}\n")

        ### VALIDATION SET ###
        distance2, indices2 = X_tree.query(center, (train_nb+val_nb))
        count_val = 0
        for index2 in indices2[train_nb:]:
            center_coords2 = tuple(centers[index2])
            img_name_val = dict_img_inverted[center_coords2]

            with open(f'{FOLDER}/val.txt', 'a') as f2:
                f2.write(f"{int(img_name_val)}\n")
            count_val += 1
        
        ### TEST SET ###
        distance3, indices3 = X_tree.query(center, (len(dict_img)))
        # train_to_add, test_to_add, val_to_add = [], [], []
        count_test = 0
        for index3 in indices3[(train_nb+val_nb):]:
            center_coords3 = tuple(centers[index3])
            img_name_test = dict_img_inverted[center_coords3]

            with open(f'{FOLDER}/test.txt', 'a') as f3:
                f3.write(f"{int(img_name_test)}\n")
            count_test += 1

        # Assert that all images have been written, once and only once, 
        # in one of the 3 folders
        assert (len(indices) + count_val + count_test) == len(dict_img)
    return

def split(DIR_IMAGES_GEOTIFF, NUMBER_OF_SPLITS, PERCENTAGE_TRAIN, PERCENTAGE_VAL, PERCENTAGE_TEST, margin=0.03):
    '''
    Function that splits data input into train-val-test sets, 
    making sure the 3 datasets do not overlap and respect the percentage asked (with tolerance margin).
    
    input: DIR_IMAGES_GEOTIFF = folder directory of geotif images, str
    NUMBER_OF_SPLITS = number of split wished, int
    PERCENTAGE_TRAIN = percentage of data to set as training data, float in [0-1]
    PERCENTAGE_VAL = percentage of data to set as validation data, float in [0-1]
    PERCENTAGE_TEST = percentage of data to set as test data, float in [0-1]
    margin = percentage of tolerance for data_split. Fix to zero to get the exact percentage, float in [0-1].
    
    output: a folder per split, containing 3 txt files (train, validation, test)
    '''
    # WE GATHER AND KEEP TRACK OF IMAGES' CENTERS
    dict_img = {}
    dict_img_inverted = {}
    center_x, center_y = [], []
    count = 0
    size_img = 0
    for image in os.listdir(DIR_IMAGES_GEOTIFF):
        if image[-4:] == ".tif":
            img_path = f"{DIR_IMAGES_GEOTIFF}/{image}"
            _, _, _, [ulx_deg, uly_deg, lrx_deg, lry_deg] = get_img_bbox(img_path)
            x = np.float128(ulx_deg + (lrx_deg-ulx_deg))
            y = np.float128(lry_deg + (lry_deg-uly_deg))
            # We determine the radius search around an image, according to the biggest size of the img
            if size_img == 0:
                size_img = np.amax([(lrx_deg-ulx_deg), (lry_deg-uly_deg)])
            # Assert that there are no image duplicates (name or content)
            # Duplicates will be ignored: Only a unique image (last image encountered) is kept.
            try:
                assert image[:-4] not in dict_img
                try:
                    # TWO ASSERTIONS ARE PASSED, THE IMAGE IS KEPT
                    assert (x,y) not in dict_img_inverted
                    center_x.append(x)
                    center_y.append(y)
                    dict_img[image[:-4]] = (x, y)
                    dict_img_inverted[(x,y)] = image[:-4]
                    count += 1
                except:
                    continue
                    print(f"All images are not unique, check out your input. E.g. {dict_img_inverted[(x,y)]}.tif and {image[:-4]}.tif are similar.")
            except:
                continue
                print(f"All images do not have a unique name, check out your input. E.g. {image[:-4]}.tif exists at least twice.")
    try:
        assert count == len(dict_img)
    except:
        print("Please, check out your input images.")

    # KD-TREE FROM IMG CENTERS
    centers = np.array(list(zip(center_x,center_y)))
    X_tree = scipy.spatial.KDTree(centers)

    # HOW MANY TRAIN/VAL/TEST SAMPLES
    train_nb = round(len(dict_img) * PERCENTAGE_TRAIN)
    val_nb = round((len(dict_img) - train_nb) * PERCENTAGE_VAL/(PERCENTAGE_VAL+PERCENTAGE_TEST))
    test_nb = (len(dict_img) - (val_nb + train_nb))
    print(f"There are {len(dict_img)} images. Split has {train_nb} train. images, {val_nb} val. images, {test_nb} test images, with {(margin*100)} percents precision.")
    assert (train_nb + val_nb + test_nb == len(dict_img))

    for i in range(NUMBER_OF_SPLITS):

        ### 1 ### CREATE DATA-SPLIT, SAVED IN SETS()

        # WE ENSURE THAT OVERLAPPING IMAGES ARE IN THE SAME DATASET USING A QUEUE STRUCTURE
        train_data = set()
        radius_search = (size_img/2) * math.sqrt(2) # diagonal of the image
        # First while loop to reach the exact number of training data +/- margin
        while len(train_data) < (train_nb - margin*train_nb) or len(train_data) > (train_nb + margin*train_nb):
            # RESET ALL VALUES TO TRY AGAIN IF TRAIN_NUMBER NOT WITHIN PERCENTAGE +/- MARGIN WISHED
            img_queue = []
            train_data = set()
            keep_track = set()
            first_train_point = random.sample(list(dict_img), 1)[0]
            img_queue.append(first_train_point)
            keep_track.add(first_train_point)
            # Second while loop to obtain training data until there are no more overlaps with other datasets
            while len(train_data) < (train_nb - margin*train_nb) or len(img_queue) > 0:
                # If there are no image in the queue, we select a random one to add in the queue
                if len(img_queue) == 0:
                    random_center = random.sample(list(dict_img), 1)[0]
                    while random_center in keep_track:
                        random_center = random.sample(list(dict_img), 1)[0]
                    img_queue.append(random_center)
                    keep_track.add(random_center)
                else:
                    point = img_queue.pop(0)
                    assert point in keep_track
                    train_data.add(point)  
                    # WE RETRIEVE IMAGE(S) THAT MIGHT OVERLAP USING THE KD-TREE WITH RADIUS SEARCH OF IMAGE DIAGONAL
                    center_x, center_y = dict_img[point]
                    indices2 = X_tree.query_ball_point(x=[center_x, center_y], r=radius_search)
                    # image1 = dict_img_inverted[(center_x, center_y)]
                    # img_path1 = f"{DIR_IMAGES_GEOTIFF}/{image1}.tif"
                    # poly1 = get_poly_from_img(img_path1)
                    for index in indices2:
                        image_to_queue = dict_img_inverted[tuple(centers[index])]
                        # img_path2 = f"{DIR_IMAGES_GEOTIFF}/{image_to_queue}.tif"
                        # poly2 = get_poly_from_img(img_path2)
                        # if image_to_queue not in keep_track and poly1.intersects(poly2):
                        if image_to_queue not in keep_track:
                            img_queue.append(image_to_queue)
                            keep_track.add(image_to_queue)
                        else:
                            continue
            if len(train_data) < int(train_nb - margin*train_nb) or len(train_data) > int(train_nb + margin*train_nb):
                print(f"The random split is run again, because we obtained {len(train_data)} training data samples instead of {int(train_nb - margin*train_nb)} to {int(train_nb + margin*train_nb)}.")

        assert(len(img_queue) == 0 and len(train_data) >= train_nb)

        # We do the same with validation data for the remaining images
        keep_track_bak = keep_track.copy()
        val_data = set()
        # First while loop to reach the exact number of validation data +/- margin
        while len(val_data) < (val_nb - margin*val_nb) or len(val_data) > (val_nb + margin*val_nb):
            # RESET ALL VALUES TO TRY AGAIN
            img_queue = []
            val_data = set()
            keep_track = keep_track_bak.copy()
            assert (len(keep_track) == len(train_data))
            first_val_point = random.sample(list(dict_img), 1)[0]
            while first_val_point in keep_track:
                first_val_point = random.sample(list(dict_img), 1)[0]
            img_queue.append(first_val_point)
            keep_track.add(first_val_point)
            # Second while loop to obtain validation data until there are no more overlaps with other datasets
            while len(val_data) < (val_nb - margin*val_nb) or len(img_queue) > 0:
                # If there are no image in the queue, we select a random one to add in the queue
                if len(img_queue) == 0:
                    random_center2 = random.sample(list(dict_img), 1)[0]
                    while random_center2 in keep_track:
                        random_center2 = random.sample(list(dict_img), 1)[0]
                    img_queue.append(random_center2)
                    keep_track.add(random_center2)
                else:
                    point2 = img_queue.pop(0)
                    assert point2 in keep_track
                    val_data.add(point2)  
                    # WE RETRIEVE IMAGE(S) THAT MIGHT OVERLAP
                    center_x, center_y = dict_img[point2]
                    indices2 = X_tree.query_ball_point(x=[center_x, center_y], r=radius_search)
                    for index2 in indices2:
                        image_to_queue2 = dict_img_inverted[tuple(centers[index2])]
                        if image_to_queue2 not in keep_track:
                            img_queue.append(image_to_queue2)
                            keep_track.add(image_to_queue2)
                        else:
                            continue
            if len(val_data) < (val_nb - margin*val_nb) or len(val_data) > (val_nb + margin*val_nb):
                print(f"The random split is run again, because we obtained {len(val_data)} training data samples instead of {int(val_nb - margin*val_nb)} to {int(val_nb + margin*val_nb)}.")

        assert(len(keep_track) == (len(train_data) + len(val_data)))
        # THE REMAINING POINTS ARE PUT INTO TEST DATASET
        test_data = set()
        for img_name in dict_img:
            if img_name not in keep_track:
                test_data.add(img_name)
        
        assert(len(img_queue) == 0) \
            and (len(train_data) + len(val_data) + len(test_data)) == len(dict_img)

        ### 2 ### WRITE FINAL DATA SPLIT TO TXT FILES

        # create folder structure, if not existing
        print(f"Split {i+1} contains {len(train_data)} train. data, {len(val_data)} val. data, {len(test_data)} test data.")
        FOLDER = DATA_SPLIT_FOLDER + '_' + str(i+1)
        if os.path.isdir(FOLDER):
            shutil.rmtree(FOLDER)
        if not os.path.isdir(FOLDER):
            os.mkdir(FOLDER)

        with open(f'{FOLDER}/train.txt', 'w') as f_train:
            for training_img in train_data:
                f_train.write(f"{int(training_img)}\n")
        with open(f'{FOLDER}/val.txt', 'w') as f_val:
            for val_img in val_data:
                f_val.write(f"{int(val_img)}\n")
        with open(f'{FOLDER}/test.txt', 'w') as f_test:
            for test_img in test_data:
                f_test.write(f"{int(test_img)}\n")

        print(f"Split files are written for split {i+1} in {FOLDER}.")

    return

if create_split == True:
    split(DIR_IMAGES_GEOTIFF, NUMBER_OF_SPLITS, \
    PERCENTAGE_TRAIN, PERCENTAGE_VAL, PERCENTAGE_TEST)

if visualise_split == True:
    datasets_to_geojson(DATA_SPLIT_FOLDER_TO_VISUALIZE, DIR_IMAGES_GEOTIFF)