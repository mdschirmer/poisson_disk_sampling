###############
# Introduction
###############

This is a software package which provides Poisson Disk Sampling. For details see 'docs'.

###############
# Install
###############

Use:

python setup.py install

Requirements which cannot be satisfied by the setup script is graph-tool (see below).

###############
# GRAPHTOOL
###############

It is highly recommended that graph-tool is compiled with its multi-core capabilities activated. This can be done, e.g. by compiling graph-tool from source, using 

./configure --enable-openmp
make -jX  # X is the number of cores used in making
make install

Graph-tools relies on a variety of libraries. These are the libraries I had to install on Ubuntu 16.04, in order for the configuration to proceed without errors. (no guarantee of completeness; date: 2016-05-02)

============
sudo apt install libboost-all-dev
pip install scipy
sudo apt install libcgal-dev
sudo apt install libcairo2-dev
pip install matplotlib
sudo apt install python-cairo-dev
sudo apt install libsparsehash-dev
============