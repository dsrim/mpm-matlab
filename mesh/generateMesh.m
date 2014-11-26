function [c4n,n4e] = generateMesh(xlim,ylim,nx,ny)

xx = linspace(xlim(1),xlim(2),nx+1); 
yy = linspace(ylim(1),ylim(2),ny+1);
[X,Y] = meshgrid(xx,yy);
nxy = reshape(1:((nx+1)*(ny+1)),nx+1,ny+1)';
c4n = [X(:),Y(:)];
nNodes = size(c4n,1);
nxy1 = nxy(1:(end-1),1:(end-1));
nxy2 = nxy(2:end,1:(end-1));
nxy3 = nxy(2:end,2:end);
nxy4 = nxy(1:(end-1),2:end);
n4e = [nxy1(:),nxy2(:),nxy3(:),nxy4(:)];
  
end