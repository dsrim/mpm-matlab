function [x4pn,e4pn,v4pn,F4pn,J4p,rho4pn,vol4pn,sigma4pn] = I2pP1(...
                            x4p,v4p,m4p,e4p,c4n,n4e,e4n,n4s,s4e,rhop0,...
                            vI,vIold,F4p,nrPts,nrNodes,dt,lambda,mu)
% update material pts/grid values
F4pn = zeros(2,2,nrPts);
J4p = zeros(nrPts,1);
x4pn = zeros(nrPts,2);
v4pn = zeros(nrPts,2);
sigma4pn = zeros(2,2,nrPts);
for p = 1:nrPts
  % fix p then find nonzero NIxp
  nodes = n4e(e4p(p,:),:);
  c4nloc = c4n(nodes,:); 
  vIloc = vI(nodes,:); 
  vIoldloc = vIold(nodes,:);
%   NIxp = shapeP1(x4p(p,:),c4nloc);
  dNIxp = shapeP1g(c4nloc);
  dNIxp = repmat(NIxp,1,2);
  F4pn(:,:,p) = (eye(2) + dt*vIloc'*dNIxp)*F4p(:,:,p);
  x4pn(p,:) = x4p(p,:) + dt*NIxp'*vIloc;
  v4pn(p,:) = v4p(p,:) + NIxp'*(vIloc-vIoldloc);
  J4p(p) = det(F4pn(:,:,p));
  if J4p(p) < 0
      display(['warning: negative Jacobian for pt: ' num2str(p)])
  end
  sigma4pn(:,:,p) = (lambda*log(J4p(p))*eye(2) ...
                       + mu*(F4pn(:,:,p)*F4pn(:,:,p)' - eye(2)))./J4p(p);
end
rho4pn = rhop0./J4p;
vol4pn = m4p./rho4pn;
e4pn = updateE4pS(c4n,n4e,e4n,x4p,e4p);
end
