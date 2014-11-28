function vI = normalBC(vI,bdNormals)
N = length(vI);
vI(bdNormals(:,1) + (abs(bdNormals(:,2))-1)*N) = 0;
end