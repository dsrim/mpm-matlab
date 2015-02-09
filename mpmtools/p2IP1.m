function [mI,mvI,mIJ,vI,fI] = p2IP1(x4p,e4p,c4n,n4e,n4s,s4e,e4s,m4p,vol4p,v4p,...
                                    b4p,sigma4p,bdNormals)
% map values on material points to the grid
nrNodes = size(c4n,1);
mI = zeros(nrNodes,1);
mIJ = sparse(zeros(nrNodes,nrNodes));
mvI = zeros(nrNodes,2);
fI = zeros(nrNodes,2);
nrPts = size(x4p,1);
for p = 1:nrPts
  % Fix p then find I to contribute to
  nodes = n4e(e4p(p,:),:);
  NIxp = shapeP1(x4p(p,:),c4n(nodes,:));
  dNIxpdx = shapeP1g(c4n(nodes,:));
%   dNIxpdx = repmat(NIxp,1,2);
  NIxNI = NIxp*NIxp';
  NIxNI = [diag(NIxNI); diag(NIxNI,1); diag(NIxNI,2)];
  ind1 = [1 2 3 1 2 1]; ind2 = [1 2 3 2 3 3];
  mIJ = mIJ + sparse(nodes(ind1),nodes(ind2),m4p(p)*NIxNI,nrNodes,nrNodes);
  col = ones(length(nodes),1); cols = [col; 2*col];
  mI = mI + sparse(nodes, col, m4p(p)*NIxp,nrNodes,1);
  mvI = mvI + sparse([nodes;nodes],cols, m4p(p).*NIxp*v4p(p,:),nrNodes,2);
  fI = fI + sparse([nodes;nodes],cols,-vol4p(p)* dNIxpdx*sigma4p(:,:,p)'...
      + m4p(p).*NIxp*b4p(p,:),nrNodes,2);
end
mvI(mI < 1e-14) = 0;
mvI = normalBC(mvI,bdNormals);
fI = normalBC(fI,bdNormals);
vI = mvI./repmat(mI,1,2);
end
