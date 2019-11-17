# Libraries
# =========

# Standars
import random, sys

# Image Processing (general)
import cv2
import numpy as np

# Image Processing (watershed segmentation)
from skimage.feature import peak_local_max
from skimage.morphology import watershed
from scipy import ndimage
from scipy.sparse import csr_matrix
from scipy.ndimage.filters import gaussian_filter

def get_cnt_mask(cluster_index, sp_arr, labels_shape):
    # time.sleep(1)
    cnt_y_index, cnt_x_index = np.unravel_index(sp_arr.indices[sp_arr.data == cluster_index], labels_shape)
    
    cnt_x_min, cnt_x_max, cnt_y_min, cnt_y_max = np.min(cnt_x_index), np.max(cnt_x_index), np.min(cnt_y_index), np.max(cnt_y_index)
    cnt_topLeft_P = (cnt_x_min, cnt_y_min)
    cnt_img_h , cnt_img_w = cnt_y_max - cnt_y_min + 1 , cnt_x_max - cnt_x_min + 1
    cnt_img_shape = (cnt_img_h, cnt_img_w)
    cnt_mask_x_index = cnt_x_index - cnt_x_min
    cnt_mask_y_index = cnt_y_index - cnt_y_min
    cnt_mask_xy_index = (cnt_mask_y_index , cnt_mask_x_index)
    cnt_mask = np.zeros(cnt_img_shape, dtype="uint8")
    cnt_mask[cnt_mask_xy_index] = 1    
    
    
    return cnt_mask, cnt_topLeft_P


def get_contour_in_mask(cnt_mask, cnt_topLeft_P):
    
    # detect contours in the cnt_mask and grab the largest one
    cnts = cv2.findContours(cnt_mask.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)[-2]
    c = max(cnts, key = cv2.contourArea)
    c =  c.reshape(c.shape[0],c.shape[2])
    
    # apply offset to account for the correct location of the cell in the image
    c[:,0] += cnt_topLeft_P[0]
    c[:,1] += cnt_topLeft_P[1]
    
    return c

def random_color():
    r = lambda: random.randint(0,255)
    return (r(),r(),r())

def fill_holes_in_binary_image(img_binary):
    
    '''
    THIS FUNCTION NEEDS TO BE OPTIMISED.
    
    ** SLOW MODE **
    
    '''

    # find contours
    contour,hier = cv2.findContours(img_binary,cv2.RETR_CCOMP,cv2.CHAIN_APPROX_SIMPLE)

    # black blank image
    img_binary_new = np.zeros(shape=img_binary.shape, dtype=np.uint8)

    # refill contours
    for cnt in contour:
            cv2.drawContours(img_binary_new,[cnt],0,255,-1)

    return img_binary_new

def downsize_image(img_8bit, downsize_factor):
    
    # create a deep copy
    img_8bit_downsized = img_8bit.copy()
    
    # down sizing
    for i in range(downsize_factor): img_8bit_downsized = cv2.pyrDown(img_8bit_downsized)
        
    # display message
    table = {'input size': str(img_8bit.shape), 'ouput size': str(img_8bit_downsized.shape)}
    print(f'Image is downsized by a factor of {downsize_factor}.')
    for name, size in table.items():
        print(f'{name:10} ==> {size:20}')
    
    return img_8bit_downsized

# Segmentation (FUNCTION)
# -----------------------

def apply_w_shed_segmentation(img_8bit, adapThresh_blcokSize, adapThresh_constant, min_distance):

    
    img_8bit_copy = img_8bit.copy()
    img_8bit_copy = cv2.cvtColor(img_8bit_copy, cv2.COLOR_GRAY2BGR)

    
    height, width = img_8bit.shape[:2]


    # binarise image
    img_binary = cv2.adaptiveThreshold(img_8bit, 255,
                                       cv2.ADAPTIVE_THRESH_MEAN_C,
                                       cv2.THRESH_BINARY,
                                       adapThresh_blcokSize,
                                       adapThresh_constant)
    
    # modify the binary image by filling holes
    img_binary = fill_holes_in_binary_image(img_binary)

    
    # watershed segmentation
    D = ndimage.distance_transform_edt(img_binary)
    localMax = peak_local_max(D, indices=False, min_distance=min_distance, labels=img_binary)
    markers = ndimage.label(localMax, structure=np.ones((3, 3)))[0]
    
    labels = watershed(-D, markers, mask=img_binary) # If returned warning, upgrade to latest Skimage    
    sp_arr = csr_matrix(labels.reshape(1,-1))
    labels_shape = labels.shape
    cluster_index_total = len(np.unique(sp_arr.data))
    cluster_index_list = np.arange(1, cluster_index_total + 1 ) # excluding index = [background]
    
    # keep positions
    points = np.zeros([20000, 2], dtype=int)
    X , Y = 0 , 1
    n_cell = 0
    
    # go through non-zero values in the sparsed version of the 'labels' array
    # for i in np.unique(sp_arr.data)[1:10]:  # as a bonus this `unique` call should be faster too
    for cluster_index in cluster_index_list: # OR for i in np.unique(sp_arr.data)

        # get mask image and its top-left corner position
        cnt_mask , cnt_topLeft_P = get_cnt_mask(cluster_index, sp_arr, labels_shape)

        # detect contours in the cnt_mask and grab the largest one
        cnt = get_contour_in_mask(cnt_mask, cnt_topLeft_P)

        if len(cnt) >= 5:
            
            ellipse = cv2.fitEllipse(cnt)
            (x, y), (MA, ma), angle = ellipse
            centre = (int(x),int(y))
            
            points[n_cell][X] = x
            points[n_cell][Y] = height - y
            
            ellipse_area = np.pi * (MA/2.0) * (ma/2.0)
#             cv2.ellipse(img_8bit_copy, ellipse, (255,0,0), 1)
            cv2.drawContours(img_8bit_copy, [cnt], 0, (0,255,0), 1)
#             cv2.circle(img_8bit_copy, centre, 1, (255,255,255), 1)
            n_cell += 1
            
    print ("Total No of cells: ", n_cell)
    points = points[0:n_cell]
    
    return points, img_8bit_copy
