function L2error = run2DMPM(n)
% Computes Material Point Method solution and returns computed L^2 error
%
% L2error = run2DMPM(n)
% 
% n: number of elements along x-direction
% L2error: computed L^2 error between true and approx solution.
%
% Donsub Rim (drim@uw.edu)
%


%% Run 2D MPM code
addpath(genpath('./'))
xlim = [0,1]; ylim = [0,1]; 
nmp4e = 4; 
dx = (xlim(2) - xlim(1))/n;

%% Initialization
lambda = 0.3;                                   % Lame const
mu = 0.4;                                       % Lame const
c0 = sqrt(mu);                                  % prob param 
t = 0;                                          % init time
T = 2/c0 ;                                      % final time
dt = 0.5*(dx);                                  % time step
A = 0.01;                                       % prob param
v0 = @(X) A*sin(pi*X);                          % initial velocity
f = @(Z,t) A/c0*cos(pi*Z).*sin(c0*pi*t) + 1;    % defd for convenience
df = @(Z,t) -A*pi/c0*sin(pi*Z).*sin(c0*pi*t);   % defd for convenience
b = @(X,t) [(c0^2 - mu)*df(X(:,1),t) ...        % body force
  + (lambda*log(abs(f(X(:,1),t).*f(X(:,2),t))) ...
  - (lambda + mu)).*df(X(:,1),t)./f(X(:,1),t).^2, ...
  (c0^2 - mu)*df(X(:,2),t) ...
  + (lambda*log(abs(f(X(:,1),t).*f(X(:,2),t))) ...
  - (lambda + mu)).*df(X(:,2),t)./f(X(:,2),t).^2];
rho0 = @(x) 0*x(:,1)+1;                         % init density = 1
sigma0 = @(x) zeros(2,2,length(x)) ;            % init stress = 0
uexact = @(X,t) A/(c0*pi)*sin(pi*X)*sin(c0*pi*t);   % exact displacement
vexact = @(X,t) A*sin(pi*X)*cos(c0*pi*t);           % exact velocity

% Grid / Material Point set-up.
nx = n; ny = n;                                 % nr of elts along x,y
[c4n,n4e,inNodes,bdNodes,bdNormals,bdTangents,inComps,vol4e] = ...
                                        generateRectMesh(xlim,ylim,nx,ny);
e4n = computeE4n(n4e);                          % compute e4n
[x4p,e4p] = generateMp(c4n,n4e,@isbody,nmp4e);  
nrNodes = size(c4n,1);
nrPts = size(x4p,1);

% Initialize variables.
b4p = b(x4p,0);                          
v4p = v0(x4p);
vol4p = zeros(nrPts,1);                  % init var vol4p
rho4p = zeros(nrPts,1);                  % init var rho4p
parfor p = 1:nrPts
  vol4p(p) = vol4e(1)/nmp4e;             % init vol
  F4p(:,:,p) = eye(2,2);                 % init deformation gradient
end
J4p = zeros(nrPts,1);                    % init var J4p
rho4p(:,1) = rho0(x4p(:,1));             % init mass density
sigma4p = sigma0(x4p);                   % init stress
m4p(:,1) = vol4p(:,1).*rho4p(:,1);       % init material point mass
uI = zeros(nrNodes,2);                   % init nodal vars
vIold = zeros(nrNodes,2);
vI = zeros(nrNodes,2);

%% MPM Loop

while (t < (T - 100*eps))
    
if t + dt > T
     t = T;
    dt = T - t;
end

% Interpolate to Grid
[mI,mvI,fI] = p2I(x4p,e4p,c4n,n4e,m4p,vol4p,v4p,...
                  b4p,sigma4p,nrPts,nrNodes,bdNormals);

% Solve on Grid (update momentum)
posI = mI > 1e-12;    
vIold(posI,:) = mvI(posI,:)./repmat(mI(posI),1,2);  
mvI(inComps) = mvI(inComps) + dt*fI(inComps);
vI(posI,:) = mvI(posI,:)./repmat(mI(posI),1,2);     
vI(~posI,:) = zeros(sum(~posI,1),2);       % remove node values close to 0
uI = uI + dt*vI;

% Update Material Pts
t = t + dt;                       % time step
b4p = b(x4p,t);                   % update body force
[x4p,e4p,v4p,F4p,J4p,rho4p,vol4p,sigma4p] = I2p(...
                                     x4p,v4p,m4p,e4p,c4n,n4e,e4n,rho4p(:,1),...
                                     vI,vIold,F4p,nrPts,nrNodes,dt,lambda,mu);

display(['time = ' num2str(t,'%1.4f') '/' num2str(T,'%1.4f')]);


plotSol(c4n,n4e,x4p,t,v4p,vexact(x4p,t),uI,uexact(c4n,t))
end


plotSol(c4n,n4e,x4p,t,v4p,vexact(x4p,t),uI,uexact(c4n,t))
[L2error,upbd] = computeL2error(uI,uexact,c4n,n4e,T);
end

function [xx, yy] = isbody(xx,yy)
% Given (xx,yy) returns only those points that are inside the body
   j = logical((xx >= 0) .* (xx <= 1) .* (yy >= 0) .* (yy <= 1));
   xx = xx(j);
   yy = yy(j);
end

