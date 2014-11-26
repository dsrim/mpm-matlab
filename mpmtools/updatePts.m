function [x4pn,e4pn,v4pn,F4pn,J4p,rho4pn,vol4pn,sigma4pn] = updatePts(...
                            x4p,v4p,m4p,mI,volp,e4p,c4n,n4e,e4n,...
                            rhop0,vI,vIold,F4p,FI,nrPts,bdNodes,E,dt,lambda,mu)
nrNodes = size(n4e,1);
% update material pts/grid values
Fpfactor = zeros(2,2,nrPts);
F4pn = zeros(2,2,nrPts);
J4p = zeros(nrPts,1);
x4pn = zeros(nrPts,2);
v4pn = zeros(nrPts,2);
sigma4pn = zeros(2,2,nrPts);
for p = 1:nrPts
  % fix p then find nonzero NIxp
  nodes = n4e(e4p(p,:),:);
  NIxp = shapeR(x4p(p,:),c4n(nodes,:));
  dNIxp = shapeRg(x4p(p,:),c4n(nodes,:));
  Fpfactor(:,:,p) = eye(2,2) + dt*vI(nodes,:)'*dNIxp;
%   Fpfactor(:,:,p) = diag(diag(Fpfactor(:,:,p)));
  F4pn(:,:,p) = Fpfactor(:,:,p)*F4p(:,:,p);
  x4pn(p,:) = x4p(p,:) + dt*NIxp'*vI(nodes,:);
  v4pn(p,:) = v4p(p,:) + NIxp'*(vI(nodes,:)-vIold(nodes,:));
  J4p(p) = det(Fpfactor(:,:,p));
  if J4p(p) < 0
    display(['warning negative det(F) for p=' num2str(p) ])
  end
  sigma4pn(:,:,p) = (lambda*log(J4p(p))*eye(2) + mu*(F4pn(:,:,p)*F4pn(:,:,p)' - eye(2)))./J4p(p);
end

rho4pn = rhop0./J4p;
vol4pn = m4p./rho4pn;
e4pn = updateE4p(c4n,n4e,e4n,x4p,e4p);

end