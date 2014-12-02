function [x4pn,e4pn,v4pn,F4pn,J4p,rho4pn,vol4pn,sigma4pn] = I2p(...
                            x4p,v4p,m4p,e4p,c4n,n4e,e4n,rhop0,...
                            vI,vIold,F4p,nrPts,nrNodes,dt,lambda,mu)
% update material pts/grid values
F4pn = zeros(2,2,nrPts);
J4p = zeros(nrPts,1);
x4pn = zeros(nrPts,2);
v4pn = zeros(nrPts,2);
sigma4pn = zeros(2,2,nrPts);
parfor p = 1:nrPts
  % fix p then find nonzero NIxp
  nodes = n4e(e4p(p,:),:);
  NIxp = shapeR(x4p(p,:),c4n(nodes,:));
  dNIxp = shapeRg(x4p(p,:),c4n(nodes,:));
  F4pn(:,:,p) = (eye(2,2) + dt*vI(nodes,:)'*dNIxp)*F4p(:,:,p);
  x4pn(p,:) = x4p(p,:) + dt*NIxp'*vI(nodes,:);
  v4pn(p,:) = v4p(p,:) + NIxp'*(vI(nodes,:)-vIold(nodes,:));
  J4p(p) = det(F4pn(:,:,p));
  sigma4pn(:,:,p) = (lambda*log(J4p(p))*eye(2) ...
      + mu*(F4pn(:,:,p)*F4pn(:,:,p)' - eye(2)))./J4p(p);
end
rho4pn = rhop0./J4p;
vol4pn = m4p./rho4pn;
e4pn = updateE4p(c4n,n4e,e4n,x4p,e4p);
end
