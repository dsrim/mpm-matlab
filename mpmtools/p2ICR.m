function [mI,mvI,fI] = p2ICR(x4p,e4p,c4n,n4e,n4s,s4e,m4p,vol4p,vp,...
                                    b4p,sigma4p,bdNormals)
% map values on material points to the grid
nrSides = size(n4s,1);
mI = zeros(nrSides,1);
mvI = zeros(nrSides,2);
fI = zeros(nrSides,2);
nrSides = size(n4s,1);
nrPts = size(x4p,1);
for p = 1:nrPts
  % Fix p then find I to contribute to
  sides = s4e(e4p(p,:),:)';
  nodes = n4e(e4p(p,:),:);
  NIxp = shapeCR(x4p(p,:),c4n(nodes,:));
  dNIxpdx = shapeCRg(c4n(nodes,:));
  col = ones(length(sides),1);
  cols = [col; 2*col];
  mI = mI + sparse(sides, col, m4p(p)*NIxp,nrSides,1);
  mvI = mvI + sparse([sides;sides],cols, m4p(p).*NIxp*vp(p,:),nrSides,2);
  fI = fI + sparse([sides;sides],cols,-vol4p(p)*dNIxpdx*sigma4p(:,:,p)'...
      +m4p(p).*NIxp*b4p(p,:),nrSides,2);
end
mvI(mI < 1e-14) = 0;
mvI = normalBC(mvI,bdNormals);
end
