function [mI,mvI,fI] = p2I(xp,e4p,c4n,n4e,m4p,vol4p,vp,...
                                    b4p,sigma4p,nrPts,nrNodes,bdNormals)
% map values on material points to the grid
mI = zeros(nrNodes,1);
mvI = zeros(nrNodes,2);
fI = zeros(nrNodes,2);
parfor p = 1:nrPts
  % Fix p then find I to contribute to
  nodes = n4e(e4p(p,:),:)';
  Nxp = shapeR(xp(p,:),c4n(nodes,:));
  dNxpdx = shapeRg(xp(p,:),c4n(nodes,:));
  col = ones(length(nodes),1);
  cols = [col; 2*col];
  mI = mI + sparse(nodes, col, m4p(p)*Nxp,nrNodes,1);
  mvI = mvI + sparse([nodes;nodes],cols, m4p(p).*Nxp*vp(p,:),nrNodes,2);
  fI = fI + sparse([nodes;nodes],cols,-vol4p(p)*dNxpdx*sigma4p(:,:,p)'...
      +m4p(p).*Nxp*b4p(p,:),nrNodes,2);
end
mvI(mI < 1e-14) = 0;
mvI = normalBC(mvI,bdNormals);
end
