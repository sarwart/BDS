function output=generate_block_image(tracks,info,atlas,N_core)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%This function is used to generate block-image from the tractogram of
%%%%%dMRI dataset. It generate block-image using stride of 1 only

%%%input
%tracks: tractogram of a dMRI dataset (set of streamlines)
%info: metadata of dMRI dataset
%atlas: intra-block atlas/parcellation 
%N_core: number of worker/processors

%%%ouputs
%conn_lookup: lookup table for intra-block connections
%block_loc: lookup table for block-image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
N=size(atlas,1);%size of the block

%%%%%%%%Look-up table for streamlines for mapping intra-block connectivity
%create look-up table for streamlines (whether a streamline traverses a
%voxel or not)
disp('Creating Look-up Tables...');
tic
lim_x=info.ImageSize(1);lim_y=info.ImageSize(2);lim_z=info.ImageSize(3);
vox=info.PixelDimensions(1);
look_up_x=zeros(length(tracks),lim_x);
look_up_y=zeros(length(tracks),lim_y);
look_up_z=zeros(length(tracks),lim_z);

 for i=1:length(tracks)

 streamline=abs(floor(tracks{i}/vox))+1;
    for j=1:length(streamline)
        if streamline(j,1)<=lim_x && streamline(j,2)<=lim_y && streamline(j,3)<=lim_z
            look_up_x(i,streamline(j,1))=1;
            look_up_y(i,streamline(j,2))=1;
            look_up_z(i,streamline(j,3))=1;
        end
    end
 
 end

toc

%creating parallel pool cluster
p = parpool(N_core);

%%%%variable initializations 
[~,lim_x]=size(look_up_x); % dMRI x-dimension
[~,lim_y]=size(look_up_y); % dMRI y-dimension
[~,lim_z]=size(look_up_z); % dMRI z-dimension

%extracting upper-traingle of intra-block connectivity matrix
nR=max(max(max(atlas)));
conn_mask=triu(ones(nR,nR));

%initialize block-image (the image dimensions are converted into linear
%indices)
x_cord=((lim_z-N+1)*(lim_y-N+1)*(lim_x-N+1)); y_cord=length(find(conn_mask));
output=spalloc(x_cord,y_cord,round(y_cord/2));


%%%%for progress update
percent=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
final=floor(x_cord.*percent);

% Initialized the queue
q = parallel.pool.DataQueue;
% After receiving new data, update_progress() will be called
afterEach(q, @update_progress);
n_completed = 0;

%%%%Create block-image (the code creates an image using stride 1)
parfor total=1:x_cord
     
    
    %initialize intra-block connectivity matrix
    conn=zeros(nR,nR);
    
    %linear index to subscripts
    [i,j,z]=ind2sub([lim_x-N+1,lim_y-N+1,lim_z-N+1],total);
    
    %extract the block at the current location
    index_i=i:i+N-1;
    index_j=j:j+N-1;
    index_z=z:z+N-1;
    
    %find streamlines in the extracted block
    [x_search,~]=find(look_up_x(:,index_i)); [y_search,~]=find(look_up_y(:,index_j)); 
    [z_search,~]=find(look_up_z(:,index_z));
    
    %if block has streamlines
    check=unique(my_intersect(my_intersect(x_search,y_search),z_search));
    if length(check)>1
        for k=1:length(check) %number of streamlines
            
            %map intra-block connectivity
            temp_index=1:N; 
            streamline=floor(tracks{check(k)}./vox);
            condition=find(streamline(:,1)>=i & streamline(:,1)<=i+N-1 & streamline(:,2)>=j & streamline(:,2)<=j+N-1 & streamline(:,3)>=z & streamline(:,3)<=z+N-1);
            
            if length(condition)>1
                snode_x=find(streamline(condition(1),1)==index_i);
                snode_y=find(streamline(condition(1),2)==index_j);
                snode_z=find(streamline(condition(1),3)==index_z);
                enode_x=find(streamline(condition(end),1)==index_i);
                enode_y=find(streamline(condition(end),2)==index_j);
                enode_z=find(streamline(condition(end),3)==index_z);
                snode=atlas(snode_x,snode_y,snode_z); enode=atlas(enode_x,enode_y,enode_z);
                %below condition to ignore connections between adjacent nodes 
                 if enode~=0 && snode~=0  && snode~=enode && conn(snode,enode)==0 ...
                         && pdist2(streamline(condition(end),:),streamline(condition(1),:))>1.5
                      conn(snode,enode)=1; conn(enode,snode)=1;
                  end
            end
        end
    end%number of tracks
    
    %extract the upper-traiangle of the connectivity matrix
    conn=conn(find(conn_mask));
    output(total,:)=sparse(conn');
    
    %update progress
    send(q, total);
    %{
    check=(total==final);
    if length(find(check)) >0
        disp(['Voxels traversed: ' num2str(percent(find(check))*100) ' %']);
    end
    %}

end % for x_cord (total number of voxels in a dMRI data)

delete(p)

catch
    disp('Error: Check input files OR Number of specified cores');
    poolobj = gcp('nocreate');
    output=[];
    if ~isempty(poolobj)
        delete(poolobj);
    end
    
end
function update_progress(~)
    n_completed = n_completed + 1;
    check=(n_completed==final);
    if length(find(check)) >0
        disp(strcat('Voxels traversed:', num2str(round(n_completed/x_cord*100)),'%'));
    end
    
end

end