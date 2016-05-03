#!/usr/bin/env python
"""
Timestamp: 2016-05
@author: MDS

Function to create some test files

"""

import os
import sys
import numpy as np
import nibabel as nib


def create_half_sphere_matrix(radius):
    
    # estimate matrix size with padding
    matrix_size = (2*radius)+21
    
    # initialise gmm and wmm matrices
    gmm = np.ones((matrix_size, matrix_size, matrix_size))
    wmm = np.zeros((matrix_size, matrix_size, matrix_size))
    
    # define sphere of radius in wmm_matrix
    idx_list = np.meshgrid(np.arange(matrix_size), np.arange(matrix_size), np.arange(matrix_size), indexing="ij")
    x_idx = np.reshape(idx_list[0], (idx_list[0].size))
    y_idx = np.reshape(idx_list[1], (idx_list[1].size))
    z_idx = np.reshape(idx_list[2], (idx_list[2].size))
    
    # define sphere centre
    centre = np.floor(matrix_size/2.)
    
    # find distances
    distances = np.sqrt(np.power(x_idx-centre,2)+np.power(y_idx-centre,2)+np.power(z_idx-centre,2))
    
    # find everything within radius
    wmm_idx = np.where(distances<=(radius))[0]
    for each in wmm_idx:
        wmm[x_idx[each],y_idx[each],z_idx[each]] = 1
    
    rm_idx = np.where(distances>=(radius+5))[0]
    for each in rm_idx:
        gmm[x_idx[each],y_idx[each],z_idx[each]] = 0

    # take out wmm from gmm
    gmm = gmm-wmm
            
    return gmm, wmm

def main(argv):
	subject = 'testfile'

	left_hemisphere = subject + '_leftcortex.nii.gz'
	right_hemisphere = subject + '_rightcortex.nii.gz'
	wmm_file = subject + '_wmm.nii.gz'

	gmm, wmm = create_half_sphere_matrix(50)
	len_x = int(gmm.shape[0]/2.)
	left_gmm = gmm.copy()
	left_gmm[0:(len_x+1),:,:] = 0

	right_gmm = gmm.copy()
	right_gmm[(len_x-1):,:,:] = 0

	affine = np.diag([1,2,3,1])
	gimg = nib.Nifti1Image(left_gmm, affine)
	gimg.to_filename(left_hemisphere)

	gimg = nib.Nifti1Image(right_gmm, affine)
	gimg.to_filename(right_hemisphere)

	wimg = nib.Nifti1Image(wmm, affine)
	wimg.to_filename(wmm_file)

if __name__ == "__main__":
    main(sys.argv)