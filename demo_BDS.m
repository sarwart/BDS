
clc; clear;

%add the scripts to current path
addpath(strcat(pwd,'\BDS'));

tic
%Load data for cortical parcellation (atlas_base) and block-image(conn_map)
block_image_file=strcat(pwd,'\data\block_image.mat');
load(block_image_file);
atlas_file=strcat(pwd,'\data\atlas_2d.mat');
load(atlas_file);

%voxel size for normalizing streamlines and mapping block-chains-to 
%streamline (voxel size coresponds to the respective dMRI data
vox=2; 

%Select BDS parameters
method='det'; %select BDS method: 'det' or 'prob'
seeding='roi'; % select seeding method: 'roi' or 'full'
stride_size=1; % value stride for generarting block-image
theta_thresh=30;%angle threshold (range: 1-90)

%calling function based on selected parameters
if strcmp(method,'det')
%det BDS
[conn,block_chains,connec_chains]=BDS_det(atlas_base,conn_map,stride_size,seeding,theta_thresh);
elseif strcmp(method,'prob')
%prob BDS
[conn,block_chains,connec_chains]=BDS_prob(atlas_base,conn_map,stride_size,seeding,theta_thresh);
end

%covert blocks to streamlines
streamlines=blocks_to_streamlines(atlas_base,block_chains,vox,'.trk');


%to save tck file
%{
tracks.data=streamlines;
tracks.total_count=num2str(length(streamlines));
tracks.count=num2str(length(streamlines));
save_file=strcat(pwd,'\data\result.tck');
write_mrtrix_tracks (tracks, save_file);
%}
%to save trk file
header_file=strcat(pwd,'\data\header_2d.mat');%specify header file for the dMRI data 
load(header_file); 
header.n_count=length(streamlines);
save_file=strcat(pwd,'\data\result.trk');
trk_write(header,streamlines,save_file);

