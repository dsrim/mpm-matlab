function [x4pn,e4pn,v4pn,F4pn,J4p,rho4pn,vol4pn,sigma4pn] = I2pCR(...
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
  sides = s4e(e4p(p,:),:);
  nodes = n4e(e4p(p,:),:);
  NIxp = shapeCRp(x4p(p,:),c4n(nodes,:));
  dNIxp = shapeCRg(c4n(nodes,:));
  F4pn(:,:,p) = (eye(2,2) + dt*vI(sides,:)'*dNIxp)*F4p(:,:,p);
  x4pn(p,:) = x4p(p,:) + dt*NIxp'*vI(sides,:);
  v4pn(p,:) = v4p(p,:) + NIxp'*(vI(sides,:)-vIold(sides,:));
  J4p(p) = det(F4pn(:,:,p));
  if J4p(p) < 0
      display(['warning: negative Jacobian for pt ' num2str(p)])
  end
  sigma4pn(:,:,p) = (lambda*log(J4p(p))*eye(2) ...
      + mu*(F4pn(:,:,p)*F4pn(:,:,p)' - eye(2)))./J4p(p);
end
rho4pn = rhop0./J4p;
vol4pn = m4p./rho4pn;
e4pn = updateE4pS(c4n,n4e,e4n,x4p,e4p);
end
