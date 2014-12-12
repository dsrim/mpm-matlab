function [c4n,n4e,inNodes,bdNodes,bdNorms,bdTans,inComps,vol4e] ...
  = generateTriUnifRectMesh(xlim,ylim,nx,ny)
% generate rectangular mesh and corresponding data structures.
xx = linspace(xlim(1),xlim(2),nx+1); 
yy = linspace(ylim(1),ylim(2),ny+1);
[X,Y] = meshgrid(xx,yy);
X = X'; Y = Y';
nxy = reshape(1:((nx+1)*(ny+1)),nx+1,ny+1)';
c4n = [X(:),Y(:)];  % coordinates for node
nrNodes = size(c4n,1);
nxy1 = nxy(1:(end-1),1:(end-1));
nxy2 = nxy(1:(end-1),2:end);
nxy3 = nxy(2:end,2:end);
nxy4 = nxy(2:end,1:(end-1));
n4e = [nxy1(:),nxy2(:),nxy4(:); nxy3(:),nxy4(:),nxy2(:)];   % nodes for elt
bdNodes = [nxy(1,:)'; nxy(:,end); nxy(end,:)'; nxy(:,1)];
inNodes = 1:nrNodes; inNodes(bdNodes) = [];
vol4e = ones(size(n4e,1),1).*((xlim(2)-xlim(1))/nx)*((ylim(2)-ylim(1))/ny)/2;
lx = ones(nx+1,1);
ly = ones(ny+1,1);
bdNorms = [bdNodes, [-2*lx; 1*ly; 2*lx; -1*ly]];
bdTans = [bdNodes, [1*lx; 2*ly; -1*lx; -2*ly]];
inComps = normalBC(logical(ones(nrNodes,2)),bdNorms);
end
