function [xy,tri,DM] = readComsol(fname)
%扫描文件获得每个单元的节点信息，以及每个域和其单元之间的关系
fp = fopen(fname, 'r');%打开文件
for i=1:18
    fgets(fp);
end
num_pts=fscanf(fp, '%d # number of mesh points\n', [1,1]);
low_pt=fscanf(fp, '%d # lowest mesh point index\n', [1,1]);
fgets(fp);
%mesh point coordinates
xy = fscanf(fp, '%lf %lf \n', [2,num_pts]);
xy = xy';
%3 vertex
for i=1:7
    fgets(fp);
end
vtx_nd_ele=fscanf(fp,'%d # number of nodes per element\n',[1,1]);
vtx_ele=fscanf(fp,'%d # number of elements\n',[1,1]);
fgets(fp);
vtx=fscanf(fp,'%d \n',[1,vtx_ele]);
vtx=vtx';

vtx_ind=fscanf(fp,'%d # number of geometric entity indices\n',[1,1]);
fgets(fp);
vtx2=fscanf(fp,'%d',[1,vtx_ind]);
%%vtx2=vtx2+1;
vtx2=vtx2';
%3 edge
for i=1:5
    fgets(fp);
end
edg_nd_ele=fscanf(fp,'%d # number of nodes per element\n',[1,1]);
edg_ele=fscanf(fp,'%d # number of elements\n',[1,1]);
fgets(fp);
edg=fscanf(fp,'%d %d \n',[2,edg_ele]);
edg=edg+1;
edg=edg';


edg_ind=fscanf(fp,'%d # number of geometric entity indices\n',[1,1]);
fgets(fp);
edg2=fscanf(fp,'%d \n',[1,edg_ind]);
edg2=edg2'+1;


%3 triangle
for i=1:5
    fgets(fp);
end
tri_nd_ele=fscanf(fp,'%d # number of nodes per element\n',[1,1]);
tri_ele=fscanf(fp,'%d # number of elements\n',[1,1]);
fgets(fp);
tri=fscanf(fp,'%d %d %d \n',[3,tri_ele]);
tri=tri'+1;


tri_ind=fscanf(fp,'%d # number of geometric entity indices\n',[1,1]);
fgets(fp);
DM=fscanf(fp,'%d',[1,tri_ind]);
DM=DM';
fclose(fp);

end


