function L2error = run2DMPM(n)

%% Run 2D MPM code
addpath(genpath('./'))
xlim = [0,1]; ylim = [0,1]; 
nmp4e = 4; 
dx = (xlim(2) - xlim(1))/n;

%% Initialization
lambda = 0.3;
mu = 0.4;
c0 = sqrt(mu);
t = 0;                                          % initial time
T = 2.0/c0 ;                                       % final time
dt = 0.01*(dx);
A = 0.01;
v0 = @(X) A*sin(pi*X);
f = @(Z,t) A/c0*cos(pi*Z).*sin(c0*pi*t) + 1;
df = @(Z,t) -A*pi/c0*sin(pi*Z).*sin(c0*pi*t);
b = @(X,t) [(c0^2 - mu)*df(X(:,1),t) ...
  + (lambda*log(abs(f(X(:,1),t).*f(X(:,2),t))) - (lambda + mu)).*df(X(:,1),t)./f(X(:,1),t).^2, ...
  (c0^2 - mu)*df(X(:,2),t) ...
  + (lambda*log(abs(f(X(:,1),t).*f(X(:,2),t))) - (lambda + mu)).*df(X(:,2),t)./f(X(:,2),t).^2];
rho0 = @(x) 0*x(:,1)+1;
sigma0 = @(x) zeros(2,2,length(x)) ;
uexact = @(X,t) A/(c0*pi)*sin(pi*X)*sin(c0*pi*t);
vexact = @(X,t) A*sin(pi*X)*cos(c0*pi*t);

% Grid / Material Point set-up.
nx = n; ny = n;
[c4n,n4e,inNodes,bdNodes,bdNormals,bdTangents,inComps,bdElts,vol4e] = ...
        generateRectMesh(xlim,ylim,nx,ny);
e4n = computeE4n(n4e);
[x4p,e4p] = generateMp(c4n,n4e,@isbody,nmp4e);
nrNodes = size(c4n,1);
nrPts = size(x4p,1);

% problem setup
vol4p = zeros(nrPts,1);                  % volume
rho4p = zeros(nrPts,1);                  % mass density
parfor p = 1:nrPts
  vol4p(p) = vol4e(1)/nmp4e;
  F4p(:,:,p) = eye(2,2);                 % deformation gradient
end
J4p = zeros(nrPts,1);
uI = zeros(nrNodes,2);                  % u at nodes
b4p = b(x4p,0);
v4p = v0(x4p);
rho4p(:,1) = rho0(x4p(:,1));
rhop0 = rho4p(:,1);
sigma4p = sigma0(x4p);
m4p(:,1) = vol4p(:,1).*rho4p(:,1);
vIold = zeros(nrNodes,2);
vI = zeros(nrNodes,2);
%% MPM Loop

while (t < T - eps)
    
if t + dt > T
    t = T;
    dt = T - t;
end
% Interpolate to Grid
[mI,mvI,fI] = p2I(x4p,e4p,c4n,n4e,m4p,vol4p,v4p,...
                                    b4p,sigma4p,F4p,nrPts,nrNodes,bdElts,bdNormals);
% Solve on Grid
posI = mI > 1e-12;
vIold(posI,:) = mvI(posI,:)./repmat(mI(posI),1,2);
mvI(inComps) = mvI(inComps) + dt*fI(inComps);
vI(posI,:) = mvI(posI,:)./repmat(mI(posI),1,2);
vI(~posI,:) = zeros(sum(~posI,1),2);
uI = uI + dt*vI;

t = t + dt;
b4p = b(x4p,t);

% Update Material Pts
[x4p,e4p,v4p,F4p,J4p,rho4p,vol4p,sigma4p] = updatePts(...
                            x4p,v4p,m4p,mI,vol4p,e4p,c4n,n4e,e4n,...
                            rhop0,vI,vIold,F4p,nrPts,nrNodes,bdElts,dt,lambda,mu);


% plotSol(c4n,n4e,x4p,t,v4p,vexact(x4p,t),vI,vexact(c4n,t))
% vIexact = vexact(c4n,t);
% % uIexact = uexact(c4n,t);
% % display(num2str(max(abs(vI(inNodes,:) - vIexact(inNodes,:)))))
display(['time = ' num2str(t,'%1.4f') '/' num2str(T,'%1.4f')]);

plotSol(c4n,n4e,x4p,t,v4p,vexact(x4p,t),uI,uexact(c4n,t))
end


[L2error,upbd] = computeL2error(uI,uexact,c4n,n4e,T);

end

function [xx, yy] = isbody(xx,yy)
   j = logical((xx >= 0) .* (xx <= 1) .* (yy >= 0) .* (yy <= 1));
   xx = xx(j);
   yy = yy(j);
end

