# **Towards Deep Learning for Connectome Mapping: A Block Decomposition Framework**

We are delighted to provide the neuroscience community with a new Block Decomposition and Stitching (BDS) framework for mapping connectomes using a diffusion-weighted Magnetic Resonance Imaging (dMRI) data. This repository provides the code that could be used to map connectomes using the two variants of the proposed BDS algorithm, namely deterministic and probabilistic versions.

The rest of this document is sectioned according to folders and scripts in this repository.

[Data](https://github.com/sarwart/BDS/tree/master/Data): The folder contains a sample 2d phantom for testing the algorithm.
- nii (2d dMRI data)
- nii (white-matter mask for ground.nii.gz)
- txt (diffusion-weight gradients and b-value used for simulating ground.nii)
- tck (5000 streamlines generated for ground.nii using [Deterministic CSD](https://mrtrix.readthedocs.io/en/latest/reference/commands/tckgen.html) with default parameters. The minimum and maximum streamline length for the generated streamlines was 10 and 50 mm respectively)
- block\_image.mat (block-image generated from ground.tck using demo\_BI.m, explained below)

[**demo\_BI.m**](https://github.com/sarwart/BDS/blob/master/demo_BI.m): provides a demonstration of generating a block-image directly from a tractogram (set of streamlines generated for a dMRI dataset) specified either in “.tck” or “.trk” format. A tractogram can be generated using any tractography algorithm. The code generates block images using a stride of 1, given that we found this stride length yielded better performance in comparison to higher strides (>1).

This script generates a block-image .mat file for a dMRI with dimension X x Y x Z (ignoring the number of gradients) using a block size of NNN with an intra-block connectivity matrix M x M. The output file is a 2d matrix with the following dimensions

((X-N+1) x (Y-N+1) x (Z-N+1), number of elements in the upper triangle of M x M)

When mapping the streamlines to block, the 3d coordinate locations of the block are converted into linear indices. The simple case of generating block-image using 2d block (NN) is demonstrated in the figure below. The same methodology is followed for converting subscripts to linear indices in MATLAB ([sub2ind](https://au.mathworks.com/help/matlab/ref/sub2ind.html?fbclid=IwAR09PjMw2w6HfAe2pVSdbW56wcrd_dQ3OmkV6DxBMDaNPZuNZ6HGFcZ35Qs))

![alt text](https://github.com/sarwart/BDS/blob/master/Data/mapping_process.png)

Auxiliary functions are stored in the [strm\_to\_block ](https://github.com/sarwart/BDS/tree/master/Tck_to_block)folder. This folder should be added to your MATLAB path.

[**demo\_BDS.m**](https://github.com/sarwart/BDS/blob/master/demo_BDS.m)[ ](https://github.com/sarwart/BDS/blob/master/demo_BDS.m)provides a demonstration of the BDS algorithm (both deterministic and probabilistic variants). It requires two input files: block-image generated using the process demonstrated in the above figure and the cortical parcellation.

Auxiliary functions are stored in the [BDS](https://github.com/sarwart/BDS/tree/master/BDS) folder. This folder should be added to your MATLAB path.
