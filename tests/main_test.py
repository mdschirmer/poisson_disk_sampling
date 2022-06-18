#!/usr/bin/env python
"""
Timestamp: 2022-06
@author: MDS

Test file for poisson disk sampling

Example:
    
    poisson_disk_sample(left_hemisphere, right_hemisphere, num_regions, wmm=wmm, output_filename='pds_results.pickle', create_niftii=True, output_basename='randomLabels-')

"""

import os
import sys
import datetime

project_dir = os.path.dirname(os.path.abspath(sys.path[0]))
sys.path.append(project_dir)
import pds

# filenames
subject = 'testfile'
left_hemisphere = subject + '_leftcortex.nii.gz'
right_hemisphere = subject + '_rightcortex.nii.gz'
wmm = subject + '_wmm.nii.gz'

# create test files (surrogate, based on sphere)
with open('tests/mk_test_files.py', "r") as source_file:
    code = compile(source_file.read(), 'mk_test_files.py', "exec")
exec(code)

# how many regions
num_regions = [100]

# parcellate and time
t0 = datetime.datetime.now()
pds.poisson_disk_sample(left_hemisphere, right_hemisphere, num_regions, wmm=wmm, output_filename='pds_results.pickle', create_niftii=True, output_basename='randomLabels-')
print(datetime.datetime.now() - t0)

