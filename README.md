# BDS
Towards Deep Learning for Connectome Mapping: A Block Decomposition Framework 

"det_BDS.m" and "prob_BDS.m" generates connectome for given block-image and nodes parcellation using Det-BDS and Prob-BDS algorithm respectively 

Example dataset provided contains DK parcellation for one dMRI ("atlas.mat" ) and corresponding block-image ("block_image.mat"). 
The block-level parcellation is given the code (4x4), which was used for mapping intra-block connectivity.
The rows of teh block-image represent the linear indices of the block's location. Columns represents the upper-triangle of symmetric intra-block connectivity matrix
Note: In the provided block-image, intra-block connections seperated by a euclidean distance of 2 are ignored 
	
