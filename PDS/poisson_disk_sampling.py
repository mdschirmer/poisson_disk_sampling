"""
Timestamp: 2016-05
@author: MD Schirmer

Function to parcellate a given volume randomly into approximately equal 
sized regions based on Poisson disk sampling.

See documentation.

Citation:
    R Bridson. Fast Poisson disk sampling in arbitrary dimensions. 
    In ACM SIGGRAPH, volume 2007, page 5, 2007.

Input:
    - left_hemisphere (grey matter mask of left hemisphere)
    - right_hemisphere (grey matter mask of right hemisphere)
    - num_regions (total number of regions to parcellate brain into)
    - wmm (optional - white matter mask, used for skeletonising grey matter)
    - output_filename (optional - file name to which the results will be saved, needs to end with pickle)
    - create_niftii (optional - flag to turn on/off the output of nii.gz parcellation files)
    - output_basename (optional - if create_niftii==True, then use this as filename base for random parcellation)

Output:
    - files:
        - spatial neighbourhood graph of voxels
        - random parcellation
        - optional: skeletonised grey matter mask, if white matter mask was given
        - left hemisphere parcellation (if create_niftii=True)
        - right hemisphere parcellation (if create_niftii=True)
        - random parcellation file (if create_niftii=True)

Example:
    
    poisson_disk_sample(left_hemisphere, right_hemisphere, num_regions, wmm=None, output_filename='pds_results.pickle', create_niftii=True, output_basename='randomLabels-')

"""

#===============================================================================
# Metadata
#===============================================================================
__author__ = 'mds'
__contact__ = 'software@markus-schirmer.com'
__date__ = '2016-05'
__version__ = '1.0'

#===============================================================================
# Import statements
#===============================================================================

import os
import nibabel as nib
import numpy as np
import graph_tool as gt
import graph_tool.topology as gtt
import pickle
import sys
import misc

#===============================================================================
# Main method to poisson disk sample a brain
#===============================================================================

def poisson_disk_sample(left_hemisphere, right_hemisphere, num_regions, wmm=None, output_filename='pds_results.pickle', create_niftii=True, output_basename='randomLabels-'):
    
    # if left and right hemisphere is given as a string (filename), load
    if ((not isinstance(left_hemisphere, str)) or (not isinstance(right_hemisphere, str))):
        raise TypeError('Input error: poisson_disk_sampling(...) expects the first two arguments to be filenames.')
    
    # check if num regions is a list
    if (not isinstance(num_regions, list)):
        raise TypeError('Input error: poisson_disk_sampling(...) expects number of regions to be a list.')
    
    # ensure output_filename ends with pickle
    if ((not (output_filename[-7::]=='.pickle'))|(not isinstance(output_filename, str))):
        raise TypeError('Input error: poisson_disksampling(...) expects output name (string) ending in ".pickle".')
    
    for which, each in enumerate(num_regions):
        if ((each%2)!=0):
            print('Input error: Number of regions uneven. Fixing it by subtracting 1...')
            num_regions[which] = each-1
        
    # set graph file names
    graph_file_left = 'graph_skel_left.xml.gz'
    graph_file_right = 'graph_skel_right.xml.gz'

    # prepare imaging data for generating voxel level graphs
    print('Preparing data...')
    imgData_left, matrix_coords_left = misc.prepare_data(left_hemisphere,wmm=wmm)
    imgData_right, matrix_coords_right = misc.prepare_data(right_hemisphere,wmm=wmm)

    # if graph files do not exist, calculate and save, else load
    if ((not os.path.isfile(graph_file_left)) & (not os.path.isfile(graph_file_right))):        
        # left hemisphere
        graph_left = generate_graph(matrix_coords_left, imgData_left.copy())
        graph_left.save(graph_file_left)
        
        #right hemisphere
        graph_right = generate_graph(matrix_coords_right, imgData_right.copy())
        graph_right.save(graph_file_right)
            
    else:
        graph_left = gt.load_graph(graph_file_left)
        graph_right = gt.load_graph(graph_file_right)
    
    print('Finding region centres...')
    if os.path.isfile(output_filename):
        with open(output_filename,'rb') as outputfile:
            [regions, thresholds_left, thresholds_right, ignore] = pickle.load(outputfile)
    else:
        regions = []
        thresholds_left = []
        thresholds_right = []
        
    # for each entry in the num_regions list
    for idx, each in enumerate(num_regions):
        print('Parcellating for %04d regions' %each)
        # find parcellations in both hemispheres
        thresh_left, regions_left = parcellate_hemisphere(graph_left, matrix_coords_left, int(each/2.))
        # right hemisphere takes in distance estimate from left hemisphere, as a starting point
        thresh_right, regions_right = parcellate_hemisphere(graph_right, matrix_coords_right, int(each/2.), thresh_left)

        # restructure to save
        regions.append([regions_left[0], regions_right[0]])
        thresholds_left.append(thresh_left)
        thresholds_right.append(thresh_right)

        print('Saving progress...')
        # save results to file
        with open(output_filename,'wb') as outputfile:
            pickle.dump([regions, thresholds_left, thresholds_right, num_regions], outputfile)
    
        if create_niftii:
            print('Creating niftii files...')
            create_parcellation(output_filename, [each], left_hemisphere, right_hemisphere, output_basename)


###############################################################################
# Poisson disk sampling functions
###############################################################################
 
def parcellate_hemisphere(graph, matrix_coords, num_regions, minR=None, max_iter = 5000):

    # estimate initial distance
    if (minR is None):
        minR = estimate_distance(num_regions, graph.num_vertices())

    num_centres = 0
    break_counter = 0
    regions_done = np.inf
    region_list = []
    while ((num_centres!=num_regions)&(break_counter<max_iter)):
        # Poisson disk sampling on left and right hemisphere to find centre points
        sample_points = find_centres(graph, minR)

        # determine number of centres
        num_centres = len(sample_points)
        #print(str(num_centres) + ' @ ' + str(minR))

        if (num_centres == num_regions):
            # write down region centres
            region_list = [matrix_coords[sample_points]]
            regions_done = minR
        else:
            # change distance: if num_centres>num_regions, then increase distance, else decrease
            if (num_centres>num_regions):
                minR = minR + 1./num_centres
            else:
                minR = minR - 1./num_centres

        # safety break for while loop. stop after 500 iterations
        break_counter = break_counter + 1

    if (break_counter==max_iter):
        print('Maximum iterations reached for %03d regions. Number of regions not found.' %num_regions)
        
    return regions_done, region_list

## Function to estimate distance between region centres (see Tammes problem, upper bound)    
def estimate_distance(num_regions, num_voxel):
    # calculate radius of sphere corresponding to a surface of size num_voxel
    r_estimate = np.power(num_voxel/(4*np.pi),1./2.)
    
    # distance estimate based on upper bound of Tammes problem (Fejes Toth,1943)
    # see documentation
#    distance = r_estimate*np.sqrt(4. - (1./np.power(np.sin(np.pi*num_regions/(6.*(num_regions-2))),2.)))
    cn = np.cos((np.pi * num_regions)/(3*num_regions -6))
    angle = np.arccos(cn/(1-cn))
    distance = np.sqrt(2*np.power(r_estimate,2)*(1-np.cos(angle)))
    
    return distance

# Poisson disk sampling as in 
# R Bridson. Fast Poisson disk sampling in arbitrary dimensions. In ACM SIGGRAPH, volume 2007, page 5, 2007.

def find_centres(graph, minR):
    edge_weights = graph.edge_properties['weights']
    num_voxel = graph.num_vertices()
    sample_points = [];
    
    # first sample point
    pointID = np.random.randint(num_voxel)
    sample_points.append(pointID)

    # list of voxels to process
    process_list = np.arange(num_voxel)

    # get all distances from sample point
    dmap = gtt.shortest_distance(graph, graph.vertex(pointID)).get_array()
    d = np.asarray(dmap)

    # mark as not to process, if distances is below minR
    process_list[d<=minR] = -1
    # create list of possible next points between minR and 2*minR
    next_pick = process_list[(d>minR) & (d<2*minR) & (process_list>0)]

    while (next_pick.size>0):

        # pick next one randomly
        randElt = np.random.randint(next_pick.size)
        pointID = next_pick[randElt]

        # add to list
        sample_points.append(pointID)

        # get distances
        dmap = gtt.shortest_distance(graph, graph.vertex(pointID), weights=edge_weights).get_array()
        d = np.asarray(dmap)

        # update list to draw next centre point from
        process_list[d<=minR] = -1
        add_to_watchlist = process_list[(d>minR) & (d<2*minR) & (process_list>0)]
        update_watchlist = next_pick[process_list[next_pick]>0]
        next_pick = np.concatenate((update_watchlist, add_to_watchlist))

    return sample_points

def create_parcellation(input_filename, num_regions, left_hemisphere, right_hemisphere, output_basename='randomLabels-'):
    # load data    
    with open(input_filename,'rb') as input_file:
        [regions, thresholds_left, thresholds_right, ignore] = pickle.load(input_file)
        
        # find number of regions that exist
        num_regions_list = []
        for each_entry in regions:
            if ((each_entry[0] is None) or (each_entry[1] is None)):
                num_regions_list.append(0)
            else:
                num_regions_list.append(len(each_entry[0])+len(each_entry[1]))
    
        # prepare imaging data for generating voxel level graphs
        imgData_left, matrix_coords_left = misc.prepare_data(left_hemisphere)
        imgData_right, matrix_coords_right = misc.prepare_data(right_hemisphere)
    
        # define graph file names
        graph_file_left = 'graph_voxel_left.xml.gz'
        graph_file_right = 'graph_voxel_right.xml.gz'
    
        # if graph files do not exist, calculate and save, else load
        if ((not os.path.isfile(graph_file_left)) & (not os.path.isfile(graph_file_right))):
            print('- Creating voxel-level graphs...')
            graph_left = generate_graph(matrix_coords_left, imgData_left.copy())
            graph_left.save(graph_file_left)
            
            graph_right = generate_graph(matrix_coords_right, imgData_right.copy())
            graph_right.save(graph_file_right)
                
        else:
            graph_left = gt.load_graph(graph_file_left)
            graph_right = gt.load_graph(graph_file_right)
    
        
        for number_of_regions in np.unique(num_regions):
            # find index for the number of regions
            idx_list = np.where(np.asarray(num_regions_list)==number_of_regions)[0]
    
            # check if this number of regions exist and assign indices to populate hemispheres
            if (len(idx_list) == 0):
                print('Number of regions not found. Try filling the gaps or choose different number of regions.')
            else:
                for count, idx in enumerate(idx_list):
                    img = nib.load(left_hemisphere)
                    print('- Populating image...')
                    ## populate image
                    file_template = output_basename + '%03d' % number_of_regions + '_%03d' % (count+1) + '.nii.gz'
                    # left hemisphere
                    left_random_label = populate_mask(graph_left, regions[idx][0], imgData_left.copy(), matrix_coords_left)
                    outimg = nib.Nifti1Image(left_random_label, header=img.get_header(), affine=img.get_affine())
                    outimg.to_filename('left_' + file_template)
                    
                    # right hemisphere
                    right_random_label = populate_mask(graph_right, regions[idx][1], imgData_right.copy(), matrix_coords_right, offset=left_random_label.max())
                    outimg = nib.Nifti1Image(right_random_label, header=img.get_header(), affine=img.get_affine())
                    outimg.to_filename('right_' + file_template)
            
                    # combined image
                    random_label = left_random_label + right_random_label
                    outimg = nib.Nifti1Image(random_label, header=img.get_header(), affine=img.get_affine())
                    outimg.to_filename(file_template)

def populate_mask(graph, sample_points, imgData, coords, offset=0):
    # find distance of all voxels to region centres
    edge_weights = graph.edge_properties['weights']
    distance_mapping = np.zeros((len(sample_points), graph.num_vertices()))
    counter = 0
    for each_centre in sample_points:
        pointID = imgData[each_centre[0],each_centre[1],each_centre[2]]-1

        # get distances
        dmap = gtt.shortest_distance(graph, graph.vertex(pointID), weights=edge_weights).get_array()
        d = np.asarray(dmap)

        distance_mapping[counter,:] = d
        counter = counter + 1

    # assign voxels to closest centre
    centres = distance_mapping.argmin(axis=0)
    for each_voxel in np.arange(len(centres)):
        imgData[coords[each_voxel][0],coords[each_voxel][1],coords[each_voxel][2]] = centres[each_voxel]+offset+1

    return imgData
    
##############################################################################
# Helper functions        
##############################################################################

# generate graph from matrix
def generate_graph(coords, imgData):
    # create graph with number of nodes = number of voxel
    num_voxel = len(coords)
    graph = gt.Graph(directed=False)
    graph.add_vertex(num_voxel)

    # create edge weight property
    edge_weights = graph.new_edge_property('double')

    # create neighbourhood indexing template and weights (spatial distance)
    all_shifts = []
    weight = []
    idx_shift = [-1, 0, 1]

    for idx_x in idx_shift:
        for idx_y in idx_shift:
            for idx_z in idx_shift:
                shift_by = np.asarray([idx_x, idx_y, idx_z])
                all_shifts.append(shift_by)
                weight.append(np.sqrt(np.sum(shift_by**2)))

    # for each voxel find neighbours and connect them according to their distance
    img_shape = imgData.shape
    for iNode in range(num_voxel):
        node_coords = coords[iNode,:]
        within_reach = []

        # which nodes are within reach?
        for each_shift in all_shifts:
            voxel_coord = node_coords - each_shift
            if ((voxel_coord[0]<img_shape[0]) & (voxel_coord[1]<img_shape[1]) &(voxel_coord[2]<img_shape[2])):
                within_reach.append(voxel_coord)

        # for each node within reach create edge with weight
        counter = 0
        for each_node in within_reach:
            target = imgData[each_node[0],each_node[1],each_node[2]]
            if ((target>0) & (target!=iNode)):
                 e = graph.add_edge(graph.vertex(iNode), graph.vertex(target-1))
                 edge_weights[e] = weight[counter]
            counter = counter + 1

    graph.edge_properties['weights'] = edge_weights
    return graph
    
