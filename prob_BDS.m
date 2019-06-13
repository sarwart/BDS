
clc; clear;
%%%%%%Loading and Initializing Parameters%%%%%%%
tic
%input file path
atlas_path=strcat('atlas.mat');
load(atlas_path);
atlas_base=round(atlas_base);
st='block_image.mat';
load(st);
seeding='ROI'; % select seeding method: 'ROI' or 'full'


%BDS parameters
stride_size=1; % value stride for proposed algorithm
theta_thresh=20;%angle threshold
conn=zeros(max(max(max(atlas_base))),max(max(max(atlas_base)))); %output of tractography

%block-level atlas
block_size=4; %block size NxNxN
atlas=create_block_atlas_3d(block_size); %create block-level atlas
nodes=max(max(max(atlas)));
index_table=find(triu(ones(nodes,nodes)));
atlas_table=zeros(nodes,3);

%block-image size
blocks_per_x=size(atlas_base,1)-block_size+1;
blocks_per_y=size(atlas_base,2)-block_size+1;
blocks_per_z=size(atlas_base,3)-block_size+1;

%extracting connections per block
time=toc;
disp(['Loading and Initializing Parameters : ' num2str(time) ' sec.']);

%%%%Creating Look-up Tables(saves computational time in every iteration)%%%%

tic
%extracting the locations of all connections
[connec_table,block_table,dir] = find(conn_map');

%location of nodes in block-level atlas 
for i=1:nodes 
    [atlas_table(i,1),atlas_table(i,2),atlas_table(i,3)]=ind2sub([block_size,block_size,block_size],find(atlas==i)); 
end


%indexes of block-level connectivity matrix into subscripts 
conn_lookup=zeros(nodes^2,2);
for i=1:nodes^2
[conn_lookup(i,1),conn_lookup(i,2)]=ind2sub([nodes,nodes],i);
end

%look-up table for blocks location 
count=0;
block_loc=zeros(blocks_per_x*blocks_per_y*blocks_per_z,3);

for i=1:blocks_per_z
 for j=1:blocks_per_y
  for k=1:blocks_per_x
    count=count+1;
    block_loc(count,1)=k; block_loc(count,2)=j; block_loc(count,3)=i;   
   end
 end
end

%creating look-up table unit for connections to unit vector mapping 
vx=zeros(nodes,nodes);
vy=zeros(nodes,nodes);
vz=zeros(nodes,nodes);
for i=1:nodes
  for j=1:nodes
      if i~=j
         v=cal_orientation_3d(i,j,atlas);
         vx(i,j)=v(1); vy(i,j)=v(2); vz(i,j)=v(3);
      else
         vx(i,j)=0; vy(i,j)=0; vz(i,j)=0;
      end
  end
end

time=toc;
disp(['Look-up Tables Creation: ' num2str(time) ' sec.']);
tic
percent=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
final=floor(length(dir).*percent);
frst=0; total=0;
for i=1:length(dir)
    
    check=(i==final);
    if length(find(check)) >0
        disp(['Connections traversed: ' num2str(percent(find(check))*100) ' %']);
    end
    %extracting indexes of block and connnectivity within a block
    conn_index=index_table(connec_table(i));
    cx=conn_lookup(conn_index,1); 
    cy=conn_lookup(conn_index,2);
    
    ix=block_loc(block_table(i),1); 
    iy=block_loc(block_table(i),2); 
    iz=block_loc(block_table(i),3);
    
    start_indx=ix; start_indy=iy; start_indz=iz; %saving location for back-tracking
    found=0; %false condition to terminate streamline
    st=block_table(i); %initialization of streamline
    start=st; stop=st; 
    first=1; % aligining vectors for forward/back-tracking
    
    %Block positions 
    start_x=((start_indx-1)*stride_size)+1; start_y=((start_indy-1)*stride_size)+1;
    start_z=((start_indz-1)*stride_size)+1;
    
    %voxels of the block
    nz = atlas_base(start_x:start_x+block_size-1, start_y:start_y+block_size-1,start_z:start_z+block_size-1);
   
  if (strcmp(seeding,'ROI') && ~isempty(find(nz,1))) || strcmp(seeding,'full')
    
    %forward tracking
    while found==0
      if first==1 %for aliging vector in two directions for forward and back-tracking
         origin_vector=[vx(cx,cy),vy(cx,cy),vz(cx,cy)];
         raw_ix=ix; raw_iy=iy; raw_iz=iz; %for finding the neighborhood block
         first=0;
      end
      %searching for neighborhood block
      raw_ix=origin_vector(1)+raw_ix;
      raw_iy=origin_vector(2)+raw_iy;
      raw_iz=origin_vector(3)+raw_iz;
        
      %location of neighborhood block
      new_ix=round(raw_ix); new_iy=round(raw_iy);
      new_iz=round(raw_iz);

      if new_ix>0 && new_ix<=blocks_per_x && new_iy>0 && new_iy<=blocks_per_y && new_iz>0 && new_iz<=blocks_per_z
        if ~(ix==new_ix && iy==new_iy && iz==new_iz)  
           %extracting location and connections in next block
            new_block=sub2ind([blocks_per_x,blocks_per_y,blocks_per_z],new_ix,new_iy,new_iz);
            if ~isempty(find(conn_map(new_block,:),1))
                temp_conn=find(conn_map(new_block,:));
                ang_comp=inf(length(temp_conn),1);
                temp_vector=zeros(length(temp_conn),3);
                
                %vector conversions of connections in next block &
                %angle calculation w.r.t current vector
                 for multi=1:length(temp_conn)
                   temp_cx=conn_lookup(index_table(temp_conn(multi)),1); 
                   temp_cy=conn_lookup(index_table(temp_conn(multi)),2);
                   temp_vector(multi,:)=[vx(temp_cx,temp_cy),vy(temp_cx,temp_cy),vz(temp_cx,temp_cy)];
                   if (temp_vector(multi,:)*origin_vector')<0
                     temp_vector(multi,:)=-temp_vector(multi,:);
                   end
                   ang_comp(multi)=cal_angle_3d(origin_vector,temp_vector(multi,:));
                 end
                        
                 %randomly selecting a vector within specified angular threshold
                  angle_index = find(ang_comp<=theta_thresh);
                 if ~isempty(angle_index)
				   angle_temp=randi(length(angle_index));
                   new_cx=conn_lookup(index_table(temp_conn(angle_index(angle_temp))),1); 
				   new_cy=conn_lookup(index_table(temp_conn(angle_index(angle_temp))),2);
                   %updating current location
				   cx=new_cx; sy=new_cy;
				   origin_vector=temp_vector(angle_index(angle_temp),:);
                   stop=new_block;
                 else
                    found=1;  %termination condition (non of the connections satisfy angle threshold)
                 end
            else
               found=1; %termination condition (next block is empty-no connections)
            end
        end
        ix=new_ix; iy=new_iy; iz=new_iz; %updating location of current block 
      else
         found=1;  %termination condition (reached end of image)
      end
    end %forward-tracking
    
   %back-tracking (same as forward tracking but with reverse current orientation)
   ix=start_indx; iy=start_indy; iz=start_indz;
   cx=conn_lookup(conn_index,1); 
   cy=conn_lookup(conn_index,2);
   first=1; found=0;

    while found==0
      if first==1
         origin_vector=[vx(cy,cx),vy(cy,cx),vz(cy,cx)];
         raw_ix=ix;
         raw_iy=iy; raw_iz=iz;
         first=0;
      end
      raw_ix=origin_vector(1)+raw_ix;
      raw_iy=origin_vector(2)+raw_iy;
      raw_iz=origin_vector(3)+raw_iz;
      new_ix=round(raw_ix); new_iy=round(raw_iy); new_iz=round(raw_iz);
         
      if new_ix>0 && new_ix<=blocks_per_x && new_iy>0 && new_iy<=blocks_per_y && new_iz>0 && new_iz<=blocks_per_z
        if ~(ix==new_ix && iy==new_iy && iz==new_iz)  
           new_block=sub2ind([blocks_per_x,blocks_per_y, blocks_per_z],new_ix,new_iy, new_iz);
           if ~isempty(find(conn_map(new_block,:),1))
                    
              temp_conn=find(conn_map(new_block,:));
              ang_comp=inf(length(temp_conn),1);
              temp_vector=zeros(length(temp_conn),3);
              for multi=1:length(temp_conn)
                        
                  temp_cx=conn_lookup(index_table(temp_conn(multi)),1); 
                  temp_cy=conn_lookup(index_table(temp_conn(multi)),2);
                  temp_vector(multi,:)=[vx(temp_cx,temp_cy),vy(temp_cx,temp_cy),vz(temp_cx,temp_cy)];%cal_orientation_3d_opt(temp_cx,temp_cy,atlas_table);
                  if (temp_vector(multi,:)*origin_vector')<0
                     temp_vector(multi,:)=-temp_vector(multi,:);
                   end
                      ang_comp(multi)=cal_angle_3d(origin_vector,temp_vector(multi,:));
              end
               angle_index = find(ang_comp<=theta_thresh);
               if ~isempty(angle_index)
				   angle_temp=randi(length(angle_index));
                   new_cx=conn_lookup(index_table(temp_conn(angle_index(angle_temp))),1); 
				   new_cy=conn_lookup(index_table(temp_conn(angle_index(angle_temp))),2);
				   cx=new_cx; sy=new_cy;
				   origin_vector=temp_vector(angle_index(angle_temp),:);
                   stop=new_block;
                else
                  found=1;
                end
           else
             found=1;
           end
        end
        ix=new_ix; iy=new_iy; iz=new_iz;
      else
      found=1;
      end
    end%while for back-track
    
    %mapping connections using end-points of block-chain
    if start~=stop
        
    %mapping blocks to voxels for "start" (one end-point)   
    [start_indx, start_indy,start_indz]=ind2sub([ blocks_per_x, blocks_per_y, blocks_per_z],start);
    start_x=((start_indx-1)*stride_size)+1; start_y=((start_indy-1)*stride_size)+1;
    start_z=((start_indz-1)*stride_size)+1;
    
    nz = atlas_base(start_x:start_x+block_size-1, start_y:start_y+block_size-1, start_z:start_z+block_size-1);
    nz1=mode(nz(nz~=0));
    
    %mapping blocks to voxels for "stop"  (other end-point)  
    [stop_indx, stop_indy, stop_indz]=ind2sub([ blocks_per_x, blocks_per_y, blocks_per_z],stop);
    stop_x=((stop_indx-1)*stride_size)+1; stop_y=((stop_indy-1)*stride_size)+1;
    stop_z=((stop_indz-1)*stride_size)+1;
    nz=atlas_base(stop_x:stop_x+block_size-1, stop_y:stop_y+block_size-1, stop_z:stop_z+block_size-1);
    nz2=mode(nz(nz~=0));
    
    
    %connectivity matrix generation
   if ~isnan(nz1) && ~isnan(nz2) && nz1~=nz2
       
       conn(nz1,nz2)= conn(nz1,nz2)+1;
       conn(nz2,nz1)= conn(nz2,nz1)+1;
   end
    clear nz1 nz2 temp1 temp2;

    end
  end%seeding methodology
    
end%total connections
time=toc;
disp(['Connectome mapping completed : ' num2str(time) ' sec.']);
clearvars -except conn

