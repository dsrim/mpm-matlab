function [mI,vI,fiI,fxI,fI] = p2I(xp,e4p,c4n,n4e,m4p,vol4p,vp,...
                                    b4p,sigma4p,Fp,nrPts,nrNodes,bdNodes,bdNormals)
% map values on material points to the grid
mI = zeros(nrNodes,1);
vI = zeros(nrNodes,2);
fiI = zeros(nrNodes,2);
fxI = zeros(nrNodes,2);
for p = 1:nrPts
  % Fix p then find I to contribute to
  nodes = n4e(e4p(p,:),:)';
  Nxp = shapeR(xp(p,:),c4n(nodes,:));
  dNxpdx = shapeRg(xp(p,:),c4n(nodes,:));
  onecol = ones(length(nodes),1);
  twocols = [onecol; 2*onecol];
  mI = mI + sparse(nodes, onecol, m4p(p).*Nxp,nrNodes,1);
  vI = vI + sparse([nodes;nodes],twocols, m4p(p).*Nxp*vp(p,:),nrNodes,2);
  fiI = fiI + sparse([nodes;nodes],twocols,-vol4p(p)*dNxpdx*sigma4p(:,:,p)',nrNodes,2);
  fxI = fxI + sparse([nodes;nodes],twocols, m4p(p).*Nxp*b4p(p,:),nrNodes,2);
end
vI = vI./repmat(mI,1,2);


% impose dirichlet boundary condition by setting components to zero.
for k = 1:length(bdNormals)
  vI(bdNormals(k,1),abs(bdNormals(k,2))) = 0;  
%   fiI(bdNormals(k,1),abs(bdNormals(k,2))) = 0;
end
fI = fiI + fxI;
end