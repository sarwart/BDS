clc;
clear; 
close all;

%add the scripts to current path
addpath(strcat(pwd,'\strm_to_block'));

%%%%%%%%%%Block specifications%%%%%%%%%%%%%%%%
%specifiy size of block
N=4;
%create block-level atlas
atlas=create_block_atlas_3d(N);
%Number of nodes
nR=max(max(max(atlas)));
%number of cores
N_core=2;

%%%%%%%%%%%dMRI signal and its corresponding track file%%%%%%%%%%%%%%%%
disp('Loading dataset...');
tic

%provide respective dMRI data or atlas for dimension extraction,
%which was also used for extracting raw dMRI blocks for training CNN
dmri_file=strcat(pwd,'\data\ground.nii.gz');
info = niftiinfo(dmri_file);
info.PixelDimensions(1)=1;

toc
%provide tck file path
track_file=strcat(pwd,'\data\ground.tck');


%for tck file
tracks = read_mrtrix_tracks(track_file);
streamlines= cell(1,length(tracks.data));
for i=1:length(tracks.data)
streamlines(i)=tracks.data(i);
end
clear tracks;

%for trk file
%{
[header, tracks]=trk_read(track_file);
streamlines= cell(1,length(tracks));
for i=1:length(tracks)
streamlines{i}=tracks(i).matrix;
end
clear tracks;
%}

%%%%%%Generate block-image from the provided tracks and dMRI data%%%%%%%%%%%% 

%generate block_image
disp('Generating Block-image...');
tic
conn_map=generate_block_image(streamlines,info,atlas,N_core);
toc
%%%%%%%%%%%%%%save block_image%%%%%%%%%%%%%%%%%%%%%%%%%
save_path=strcat(pwd,'\data\block_image.mat');
save(save_path,'conn_map');

