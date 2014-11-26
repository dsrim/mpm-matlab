function L2error = run2DMPM
close all

%% Run 2D MPM code
addpath(genpath('./'))
xlim = [0,1]; ylim = [0,1]; 
nmpe = 4; n = 4;
dx = (xlim(2) - xlim(1))/n;

%% Initialization
lambda = 0.3;
mu = 0.4;
c0 = 6;
t = 0;                                          % initial time
T = 2/c0;                                       % final time
N = floor(n^2);
dt = T/N;
E = 1;
A = 0.1;
u0 = @(X) 0.*X(:,1);
v0 = @(X) A.*[sin(pi*X(:,1)),sin(pi*X(:,2))];
f = @(Z,t) A/c0*cos(pi*Z).*sin(c0*pi*t) + 1;
df = @(Z,t) -A*pi/c0*sin(pi*Z).*sin(c0*pi*t);
b0 = @(X,t) [(c0^2 - mu)*df(X(:,1),t) ...
  + (lambda*log(abs(f(X(:,1),t).*f(X(:,2),t))) - (lambda + mu)).*df(X(:,1),t)./f(X(:,1),t).^2, ...
  (c0^2 - mu)*df(X(:,2),t) ...
  + (lambda*log(abs(f(X(:,1),t).*f(X(:,2),t))) - (lambda + mu)).*df(X(:,2),t)./f(X(:,2),t).^2];
rho0 = @(x) 0*x(:,1)+1;
sigma0 = @(x) zeros(2,2,length(x)) ;

uexact = @(X,t) A/(c0*pi)*sin(pi*X)*sin(c0*pi*t);
vexact = @(X,t) A*sin(pi*X)*cos(c0*pi*t);

% Grid / Material Point set-up.
nx = n; ny = n;
[c4n,n4e,inNodes,bdNodes,bdNormals,bdTangents,inComps,vol4e] = ...
        generateRectMesh(xlim,ylim,nx,ny);
e4n = computeE4n(n4e);
[x4p,e4p] = generateMp(c4n,n4e,@isbody,nmpe);
x4p0 = x4p;
nrNodes = size(c4n,1);
nrPts = size(x4p,1);

% problem setup
vol4p = zeros(nrPts,1);                  % volume
rho4p = zeros(nrPts,1);                  % mass density
for p = 1:nrPts
  vol4p(p) = vol4e(1)/nmpe;
  F4p(:,:,p) = eye(2,2);                 % deformation gradient
  J4p(p) = det(F4p(:,:,p));
end
FI = zeros(nrNodes,1);                  % F at nodes
uI = zeros(nrNodes,2);                  % u at nodes
b4p = b0(x4p,0);
v4p = v0(x4p);
rho4p(:,1) = rho0(x4p(:,1));
rhop0 = rho4p(:,1);
sigma4p = sigma0(x4p);
m4p(:,1) = vol4p(:,1).*rho4p(:,1);

%% MPM Loop

for j = 1:(N+1)
% Interpolate to Grid
[mI,vI,fiI,fxI,fI] = p2I(x4p,e4p,c4n,n4e,m4p,vol4p,v4p,...
                                    b4p,sigma4p,F4p,nrPts,nrNodes,bdNodes,bdNormals);
% Solve on Grid
vIexact = vexact(c4n,t);
display(['-- ' num2str(max(abs(vI(inNodes,:) - vIexact(inNodes,:))))])
vIold = vI;
vI = vI + dt*fI./repmat(mI,1,2);

for k = 1:length(bdNormals)
  vI(bdNormals(k,1),abs(bdNormals(k,2))) = 0;  
%    fI(bdNormals(k,1),abs(bdNormals(k,2))) = 0;
end
uI = uI + dt*vI;

% Update Material Pts
[x4p,e4p,v4p,F4p,J4p,rho4p,vol4p,sigma4p] = updatePts(...
                            x4p,v4p,m4p,mI,vol4p,e4p,c4n,n4e,e4n,...
                            rhop0,vI,vIold,F4p,FI,nrPts,nrNodes,E,dt,lambda,mu);
t = t + dt;
b4p = b0(x4p,t);

% plotMesh(c4n,n4e,x4p,t)

% plotSol(c4n,n4e,x4p,t,v4p,vexact(x4p,t),vI,vexact(c4n,t))
plotSol(c4n,n4e,x4p,t,v4p,vexact(x4p,t),uI,uexact(c4n,t))
vIexact = vexact(c4n,t);
% uIexact = uexact(c4n,t);
% display(num2str(max(abs(vI(inNodes,:) - vIexact(inNodes,:)))))
display(['-- ' num2str(max(abs(vI(inNodes,:) - vIexact(inNodes,:))))])
display(['time = ' num2str(t,'%1.4f')]);

end

L2error = computeL2error(uI,uexact,c4n,n4e,T);
% Output step
display(L2error)


end

function [xx, yy] = isbody(xx,yy)
   j = logical((xx >= 0) .* (xx <= 1) .* (yy >= 0) .* (yy <= 1));
   xx = xx(j);
   yy = yy(j);
end

