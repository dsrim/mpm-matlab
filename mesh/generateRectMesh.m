function [c4n,n4e,inNodes,bdNodes,bdNorms,bdTans,inComps,vol4e] ...
  = generateRectMesh(xlim,ylim,nx,ny)
% generate rectangular mesh and corresponding data structures.
xx = linspace(xlim(1),xlim(2),nx+1); 
yy = linspace(ylim(1),ylim(2),ny+1);
[X,Y] = meshgrid(xx,yy);
X = X'; Y = Y';
nxy = reshape(1:((nx+1)*(ny+1)),nx+1,ny+1)';
c4n = [X(:),Y(:)];
nrNodes = size(c4n,1);
nxy1 = nxy(1:(end-1),1:(end-1));
nxy2 = nxy(1:(end-1),2:end);
nxy3 = nxy(2:end,2:end);
nxy4 = nxy(2:end,1:(end-1));
n4e = [nxy1(:),nxy2(:),nxy3(:),nxy4(:)];
bdNodes = [nxy(1,:)'; nxy(:,end); nxy(end,:)'; nxy(:,1)];
inNodes = 1:nrNodes; inNodes(bdNodes) = [];
vol4e = ones(size(n4e,1),1).*((xlim(2)-xlim(1))/nx)*((ylim(2)-ylim(1))/ny);
l = ones(nx+1,1);
bdNorms = [bdNodes, [-2*l; 1*l; 2*l; -1*l]];
bdTans = [bdNodes, [1*l; 2*l; -1*l; -2*l]];
inComps = normalBC(logical(ones(nrNodes,2)),bdNorms);
end