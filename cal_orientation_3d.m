function vector=cal_orientation_3d(c1,c2,atlas)

N=length(atlas);
temp1=find(atlas==c1); temp2=find(atlas==c2);

[v1(1), v1(2), v1(3)] = ind2sub([N,N,N],temp1(1));
[v2(1), v2(2), v2(3)] = ind2sub([N,N,N],temp2(1));

vector=(v1-v2)/norm(v1-v2);
end