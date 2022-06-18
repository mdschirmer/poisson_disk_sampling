###############
# Introduction
###############

This is a software package which provides Poisson Disk Sampling. For details see 'docs'.

If you find the software useful and used it in one of your publications, please cite the following:

Schirmer, M.D., 2015. Developing brain connectivity: effects of parcellation scale on network analysis in neonates (Doctoral dissertation, King's College London).

###############
# Install
###############

Use:

python setup.py install

Requirements which cannot be satisfied by the setup script is graph-tool (see below).

###############
# Usage
###############

A sample script can be found here: tests/main_test.py

In-python example:
    
    poisson_disk_sample(left_hemisphere, right_hemisphere, num_regions, wmm=None, output_filename='pds_results.pickle', create_niftii=True, output_basename='randomLabels-')

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

###############
# GRAPHTOOL
###############

It is highly recommended that graph-tool is compiled with its multi-core capabilities activated. This can be done, e.g. by compiling graph-tool from source, using 

./configure --enable-openmp

make -jX  # X is the number of cores used in making

make install

Last time I executed this, I got an error message, which was solved by running

./configure --enable-openmp CXXFLAGS='-std=c++14' 

instead.

Graph-tools relies on a variety of libraries. These are the libraries I had to install on Ubuntu 16.04, in order for the configuration to proceed without errors. (no guarantee of completeness; date: 2017-04-17)

sudo apt install libboost-all-dev libcgal-dev libcairo2-dev python-cairo-dev libsparsehash-dev cairomm-1.0

pip install scipy matplotlib
