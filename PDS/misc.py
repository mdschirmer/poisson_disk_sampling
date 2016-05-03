"""
:Summary:
    Variety of helper functions to use the other packages in the data analysis

:Description:
    
    cleanup_grey_matter_mask(input_file, output_file):
        cleans up grey matter mask by only keeping the largest connected volume, i.e. it removes small disconnected voxels

:TODO:
    finish this description
    
"""
#===============================================================================
# Metadata
#===============================================================================
__author__ = 'mds'
__contact__ = 'software@markus-schirmer.com'
__date__ = '2015-06'
__version__ = '0.1'

#===============================================================================
# Import statements
#===============================================================================
import nibabel as nib
import scipy
import numpy as np
import os
import warnings
import pickle
import struct
import array
import graph_tool as gt
import graph_tool.centrality as gtc
import graph_tool.clustering as gtclust

#===============================================================================
# Helper functions
#===============================================================================
def cleanup_grey_matter_mask(input_file, output_file):
    # load image
    img = nib.load(input_file)
    # label connected components
    labelled_img = scipy.ndimage.label(img.get_data())[0]
    # find labels
    labels = np.unique(labelled_img)
    # eliminate label 0 (no intensity)
    labels = labels[labels!=0]
    # find size of each connected component
    siz = []
    for label in labels:
        siz.append(np.sum(labelled_img==label))
    
    # only keep the largest one and set the rest to 0
    for label in labels[siz!=np.max(siz)]:
        labelled_img[labelled_img==label] = 0
    labelled_img[labelled_img==5] = 1
    
    # save file
    outfile = nib.Nifti1Image(labelled_img, img.get_header(), img.get_affine())
    outfile.to_filename(output_file)

def read_data_file(data_file, csv_format_flag=True):
    # function to read tractography data file
    conn_matrix = []
    
    try:
        if not csv_format_flag:
            with open(data_file,'rb') as conn_file:
                (r,c) = struct.unpack('ii',conn_file.read(8))
                numbers = array.array('d')
                numbers.fromfile(conn_file, r*c)
                conn_matrix = np.reshape(numbers, (r,c))
        else:
            with open(data_file,'r') as conn_file:
                for each_line in conn_file:
                    conn_matrix.append(map(float, each_line[0:-1].split(',')))
                np.asarray(conn_matrix)

        return conn_matrix

    except IOError as ioErr:
        print('File error: ' + str(ioErr))

def get_connectivity_matrix(data_folder, parcel, just_num_parcels=False):
    # function to read in the connectivity matrix from various probtrackx output
    
    # determine file format and set data folder to connectivity profile
    csv_format_flag = False
    if os.path.isdir(data_folder + '/CONNECTIVITY'):
        csv_format_flag = True
        folder_name = data_folder + '/CONNECTIVITY'
    else:
        folder_name = data_folder + '/TRACTOGRAPHY/whole_brain'

    # find number of parcels
    folder_list = sorted(os.listdir(folder_name))
    num_parcels = len(folder_list)
    if just_num_parcels:
        return num_parcels
    
    if parcel > num_parcels:
        warnings.warn('Number of parcels smaller than requested parcel number. Parcel number set to 1.')
        parcel = 0

    folder_name = folder_name + '/' + folder_list[parcel] + '/'
    if not csv_format_flag:
        folder_name = folder_name + 'tracts/anisotropy/'    

    # read data
    file_list = os.listdir(folder_name)
    conn_matrix = []
    
    if csv_format_flag:
        if len(file_list)>1:
            warnings.warn('More than one file found with raw data format. The first file will be used.')

        file_name = file_list[0]

        conn_matrix = read_data_file(folder_name + file_name, csv_format_flag)
        
    else:
        if len(file_list)==1:
            file_name = file_list[0]
            if file_name[-3:] == 'raw':
                conn_matrix = read_data_file(folder_name + file_name, csv_format_flag).transpose()
            else:
                raise Exception('Only one file, but neither raw nor csv data format. (%s)' % file_name)
        else:
            file_list.sort()
            file_name = 'manyfiles'
            for each_file in file_list:
                current_line = []
                try:
                    with open(folder_name + each_file) as conn_file:
                        for each_line in conn_file:
                            current_line.append(float(each_line[:-1].split(' ')[2]))

                except IOError as ioErr:
                    print('File error: ' + str(ioErr))

                conn_matrix.append(current_line)

            conn_matrix = np.asarray(conn_matrix)

    return conn_matrix, num_parcels, folder_list[parcel]

# Functions to skeletonise grey matter masks
def skeletonise(gmm, wmm_filename):
    # gmm - grey matter mask filename
    # wmm - white matter mask filename
    # output sGMM

    # load white matter mask
    if isinstance(wmm_filename, str):
        wmm = np.squeeze(nib.load(wmm_filename).get_data())
    else:
        wmm = wmm_filename.copy()

    # prepare output data
    gm_shape = gmm.shape
    sGMM = np.zeros(gm_shape)

    # define growth template
    all_shifts = growth_template()

    # find image centre
    start_points = np.where(wmm)

    # find a first set of white matter voxels
    pick = np.random.randint(len(start_points[0]), size=100)
    neighbourhood = []
    for each_pick in pick:
        neighbourhood.append(np.asarray([start_points[0][each_pick], start_points[1][each_pick],start_points[2][each_pick]]))

    # dilute wmm to ensure it overlaps with gmm
    wmm = dilute_mask(wmm)
    
    while neighbourhood:
        neighbour = neighbourhood.pop()
        
        for each_shift in all_shifts:
            test_voxel = neighbour+each_shift
            # ensure test_voxel is within image bounds
            if ((test_voxel[0]<gm_shape[0]) & (test_voxel[1]<gm_shape[1]) & (test_voxel[2]<gm_shape[2])):
                # if test voxel is grey matter, add to skeletonised gmm
                if gmm[test_voxel[0],test_voxel[1],test_voxel[2]]:
                    sGMM[test_voxel[0],test_voxel[1],test_voxel[2]] = 1
                # else add white matter voxel to neighbourhood
                elif wmm[test_voxel[0],test_voxel[1],test_voxel[2]]:
                    neighbourhood.append(test_voxel)
                    wmm[test_voxel[0],test_voxel[1],test_voxel[2]] = 0

    return sGMM    

def dilute_mask(mask, siz=5):
    # input: binary image mask
    # output: diluted mask (diluted by size siz)
    
    # find voxels in mask
    all_mask_voxels = np.where(mask)
    all_mask_voxels = [all_mask_voxels[0], all_mask_voxels[1], all_mask_voxels[2]]

    # dilute all voxels
    for each_coord in np.arange(len(all_mask_voxels[0])):
        # find all indices to grow
        x = all_mask_voxels[0][each_coord] + np.arange(-(siz-1)/2,(siz-1)/2+1)
        y = all_mask_voxels[1][each_coord] + np.arange(-(siz-1)/2,(siz-1)/2+1)
        z = all_mask_voxels[2][each_coord] + np.arange(-(siz-1)/2,(siz-1)/2+1)

        # assure position is within image
        x = x[((x>=0) & (x<mask.shape[0]))]
        y = y[((y>=0) & (y<mask.shape[1]))]
        z = z[((z>=0) & (z<mask.shape[2]))]
        
        mask[x[0]:x[-1], y[0]:y[-1], z[0]:z[-1]] = 1

    return mask

# prepare imaging data for parcellations
def prepare_data(gmm, wmm=None):
    # load gmm (grey matter mask)
    img_data = np.int64(nib.load(gmm).get_data())
    
    # if white matter mask is given skeletonise gmm based on wmm
    if wmm:
        img_data = skeletonise(img_data.copy(), wmm)
        tmp_nii = nib.load(gmm)
        skel_file = nib.Nifti1Image(img_data.astype(int), header=tmp_nii.get_header(), affine=tmp_nii.get_affine())
        skel_file.to_filename('skel_' + gmm)
    
    
    # label individual voxels
    original_shape = img_data.shape
    img_data = img_data.reshape((img_data.size,1))
    idx = np.where(img_data)[0]
    for element in range(len(idx)):
        img_data[idx[element]] = element+1
    img_data = img_data.reshape(original_shape[0:3])
    
    # get voxel coordinates in image space
    grid = np.mgrid[0:img_data.shape[0],0:img_data.shape[1],0:img_data.shape[2]] #x,y,z coordinates
    matrix_coords = np.transpose(np.asarray([grid[0][img_data!=0], grid[1][img_data!=0], grid[2][img_data!=0]]))
    
    return img_data, matrix_coords

# define growth template to dilute (6 neighbourhood)
def growth_template():
    all_shifts = []
    
    all_shifts.append(np.asarray([0,0,-1]))
    all_shifts.append(np.asarray([0,0,1]))
    all_shifts.append(np.asarray([-1,0,0]))
    all_shifts.append(np.asarray([1,0,0]))
    all_shifts.append(np.asarray([0,-1,0]))
    all_shifts.append(np.asarray([0,1,0]))

    return all_shifts

def get_spatial_properties(filename, hemisphere_label):
    # get the spatial properties of the parcellation file filename, add label to specify which hemisphere, e.g. 'l' or 'r'
    # get image data
    if isinstance(filename,str):
        img_data = np.int64(nib.load(filename).get_data())
    else:
        img_data = filename.copy()
    img_shape = img_data.shape

    # find the labels and define the offset
    labels = np.unique(img_data[img_data!=0].astype(int))
    label_offset = np.min(labels)

    # define neighbourhood template in which to look for each voxel
    neighbourhood = growth_template()

    # initialise output variables
    conn_matrix = scipy.sparse.lil_matrix(np.zeros((len(labels), len(labels))))
    region = []

    # loop through each label
    for each_label in labels:
        binary_image = img_data == each_label
        # find label positions
        coords = np.where(binary_image)
        region.append({})
        # define size and hemisphere of label
        region[each_label-label_offset]['size'] = len(coords[0])
        region[each_label-label_offset]['hemisphere'] = hemisphere_label
        region[each_label-label_offset]['coordinates'] = np.mean(coords,axis=1)

        # loop through each voxel with the label and add the neighbours to the connectivity matrix
        for each_voxel in np.arange(len(coords[0])):
            for each_neighbour in neighbourhood:
                if (((coords[0][each_voxel]-each_neighbour[0])<img_shape[0]) & ((coords[1][each_voxel]-each_neighbour[1])<img_shape[1]) & ((coords[2][each_voxel]-each_neighbour[2])<img_shape[2])):
                    test_voxel = img_data[coords[0][each_voxel]-each_neighbour[0],coords[1][each_voxel]-each_neighbour[1],coords[2][each_voxel]-each_neighbour[2]]
                    if ((test_voxel!=each_label) & (test_voxel!=0) ):
                        # make sure the matrix is symmetric
                        conn_matrix[test_voxel-label_offset, each_label-label_offset] = 1
                        conn_matrix[each_label-label_offset, test_voxel-label_offset] = 1
   
    return conn_matrix, region
    
def preprocess_tractography_data(matrix, left_hemisphere_file, right_hemisphere_file, output_folder='.', base_name='clean_'):
    # creates network matrix based on tractgraphy and a spatial neighbourhood matrix
    # takes in results file from tractography data and turns it into a .pickle and .mat file    
    
    conn_matrix = matrix.copy()
    conn_matrix = (conn_matrix + conn_matrix.transpose())/2.

    # find left and right hemisphere spatial adjacency graph
    l_conn_matrix, l_region = get_spatial_properties(left_hemisphere_file, 'l')
    r_conn_matrix, r_region = get_spatial_properties(right_hemisphere_file, 'r')

    # combine both hemispheres
    s_matrix = scipy.sparse.lil_matrix(scipy.sparse.bmat([[l_conn_matrix, None], [None, r_conn_matrix]])).todense()

    regions = l_region + r_region
    
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    parcellation_name = left_hemisphere_file.split('-')[1].split('.')[0]
    scipy.io.savemat(output_folder + '/' + base_name + parcellation_name + '.mat', mdict={'matrix':matrix, 's_matrix':s_matrix, 'regions':regions})
    with open(output_folder + '/' + base_name + parcellation_name + '.pickle', 'wb') as datafile:
        pickle.dump([matrix, s_matrix, regions], datafile)

    return conn_matrix, s_matrix, regions


# function to create a graph from a matrix. Note: matrix should be undirected
def create_graph_from_matrix(conn_matrix):
    # initialise undirected graph    
    g = gt.Graph(directed=False)
    # add vertices corresponding to the size of conn_matrix
    dummy = g.add_vertex(conn_matrix.shape[0])
    # find edges
    edge_list = scipy.sparse.find(scipy.sparse.triu(conn_matrix))
    # add edges to graph
    for each_edge in np.arange(len(edge_list[0])):
        dummy = g.add_edge(edge_list[0][each_edge], edge_list[1][each_edge])
    
    return g

# calculate network density
def get_network_density(conn_matrix):
    return np.sum(conn_matrix!=0)/float(conn_matrix.shape[0] * (conn_matrix.shape[0]-1))

# function to calculate measures
def calculate_measures(g, tmp_measures=None, measure_list=['BC', 'T', 'E']):
    if tmp_measures is None:    
        tmp_measures = dict((k, []) for k in measure_list)

    tmp_measures['BC'].append(np.mean(gtc.betweenness(g)[0].get_array()))
    tmp_measures['T'].append(gtclust.global_clustering(g)[0])
    tmp_measures['E'].append(np.mean(gtc.closeness(g,harmonic=True).get_array()))
    
    return tmp_measures

# function to calculate the mean and the confidence interval of an array
def mean_and_confidence_interval(data, axis=0, confidence=0.95):
    a = np.array(data).astype(float)
    m = np.mean(a, axis=axis)
    lower = np.percentile(a, confidence, axis=axis)
    upper = np.percentile(a, 100-confidence, axis=axis)
    return m, lower, upper