function L2error = run2DCRMPM
% Computes Material Point Method solution and returns computed L^2 error
%
% L2error = run2DMPM(n)
% 
% n: number of elements along x-direction
% L2error: computed L^2 error between true and approx solution.
%
% Donsub Rim (drim@uw.edu)
%

n = 6;
scale = 40;

%% Run 2D MPM code
addpath(genpath('./'))
xlim = [0,1]; ylim = [0,1]; 
nmp4e = 3; 
dx = (xlim(2) - xlim(1))/n;

%% Initialization
lambda = 10;                                   % Lame const
mu = 0.4;                                       % Lame const
c0 = sqrt(mu);                                  % prob param 
% c0 = 1;
t = 0;                                          % init time
T = 2/c0 ;                                      % final time
% T = 0.1;
dt = 0.01*(dx);                                  % time step
A = 0.001;                                       % prob param
v0 = @(X) A*sin(pi*X);                          % initial velocity
f = @(Z,t) A/c0*cos(pi*Z).*sin(c0*pi*t) + 1;    % defd for convenience
df = @(Z,t) -A*pi/c0*sin(pi*Z).*sin(c0*pi*t);   % defd for convenience
b = @(X,t) [(c0^2 - mu)*df(X(:,1),t) ...        % body force
  + (lambda*log(abs(f(X(:,1),t).*f(X(:,2),t))) ...
  - (lambda + mu)).*df(X(:,1),t)./f(X(:,1),t).^2, ...
  (c0^2 - mu)*df(X(:,2),t) ...
  + (lambda*log(abs(f(X(:,1),t).*f(X(:,2),t))) ...
  - (lambda + mu)).*df(X(:,2),t)./f(X(:,2),t).^2];
% b = @(X,t) 0*X;
rho0 = @(x) 0*x(:,1)+1;                         % init density = 1
sigma0 = @(x) zeros(2,2,length(x)) ;            % init stress = 0
uexact = @(X,t) A/(c0*pi)*sin(pi*X)*sin(c0*pi*t);   % exact displacement
vexact = @(X,t) A*sin(pi*X)*cos(c0*pi*t);           % exact velocity

% Grid / Material Point set-up.
nx = n; ny = n;                                 % nr of elts along x,y
[c4n,n4e,inNodes,bdNodes,bdNormals,bdTangents,inComps,vol4e] = ...
                                        generateTriUnifRectMesh(xlim,ylim,nx,ny);
e4n = computeE4n(n4e);  % compute e4n
s4e = computeS4e(n4e);
s4e = s4e(:,[2 3 1]);
e4s = computeE4s(n4e);
n4s = computeN4s(n4e);
s4n = computeS4n(n4e);
mid4s = computeMid4s(c4n,n4s);


[x4p,e4p] = generateMpTri(mid4s,s4e,@(xx,yy)isbody(xx,yy));
nrNodes = size(c4n,1);
nrSides = size(n4s,1);
nrPts = size(x4p,1);

bdNormal4S = bdNormals4s(bdNormals,s4n,s4e,e4s);
inComps4S = normalBC(logical(ones(nrSides,2)),bdNormal4S);


% Initialize variables.
b4p = b(x4p,0);                          
v4p = v0(x4p);
vol4p = zeros(nrPts,1);                  % init var vol4p
rho4p = zeros(nrPts,1);                  % init var rho4p
for p = 1:nrPts
  vol4p(p) = vol4e(1)/nmp4e;             % init vol
  F4p(:,:,p) = eye(2,2);                 % init deformation gradient
end
% J4p = zeros(nrPts,1);                    % init var J4p
rho4p(:,1) = rho0(x4p(:,1));             % init mass density
sigma4p = sigma0(x4p);                   % init stress
m4p(:,1) = vol4p(:,1).*rho4p(:,1);       % init material point mass
uI = zeros(nrSides,2);                   % init nodal vars
vIold = zeros(nrSides,2);
vI = zeros(nrSides,2);

%% MPM Loop

while (t < (T - 100*eps))
    
if t + dt > T
     t = T;
    dt = T - t;
end

% Interpolate to Grid
[mI,mIJ,vI,fI] = p2ICR(x4p,e4p,c4n,n4e,n4s,s4e,e4s,m4p,vol4p,v4p,...
                                    b4p,sigma4p,bdNormal4S);

% Solve on Grid (update momentum)
posI = mI > 1e-14;    
% vIold(posI,:) = mvI(posI,:)./repmat(mI(posI),1,2);  
% mvI(inComps4S) = mvI(inComps4S) + dt*fI(inComps4S);
% mvI = mvI + dt*fI;
% vI(posI,:) = mvI(posI,:)./repmat(mI(posI),1,2);     
vIold = vI;
dof1 = inComps4S(:,1);
vI(dof1,1) = vI(dof1,1) + dt*(mIJ(dof1,dof1) \ fI(dof1,1));
dof2 = inComps4S(:,2);
vI(dof2,2) = vI(dof2,2) + dt*(mIJ(dof2,dof2) \ fI(dof2,2));
vI(~posI,:) = zeros(sum(~posI,1),2);       % remove node values close to 0
uI = uI + dt*vI;

% Update Material Pts
t = t + dt;                       % time step
b4p = b(x4p,t);                   % update body force
[x4p,e4p,v4p,F4p,J4p,rho4p,vol4p,sigma4p] = I2pCR(...
                         x4p,v4p,m4p,e4p,c4n,n4e,e4n,n4s,s4e,rho4p(:,1),...
                         vI,vIold,F4p,nrPts,nrNodes,dt,lambda,mu);
                     
display(['time = ' num2str(t,'%1.4f') '/' num2str(T,'%1.4f')]);
% spy(e4p)

% plotSolExact(c4n,n4e,x4p,t,v4p,vexact(x4p,t),uI,uexact(c4n,t))
figure(1); clf;
% plotTriangulation(c4n,n4e)
hold on;
plot(x4p(:,1),x4p(:,2),'r.')
% quiver(x4p(:,1),x4p(:,2),scale*v4p(:,1),scale*v4p(:,2),0,'r');
vIexact = vexact(mid4s,t);
quiver(mid4s(:,1),mid4s(:,2),scale*vIexact(:,1),scale*vIexact(:,2),0,'g');
quiver(mid4s(:,1),mid4s(:,2),scale*vI(:,1),scale*vI(:,2),0,'b');
hold off;
end

L2error = sqrt(sum(error4eCRL2(c4n, n4e, @(x)vexact(x(:,1),t), vI(:,1))) ...
          + sum(error4eCRL2(c4n, n4e, @(x) vexact(x(:,2),t), vI(:,2))));


% plotSolExact(c4n,n4e,x4p,t,v4p,vexact(x4p,t),uI,uexact(c4n,t))
% [L2error,upbd] = computeL2error(uI,uexact,c4n,n4e,T);
% L2error = 0;
end

function [xx, yy] = isbody(xx,yy)
% Given (xx,yy) returns only those points that are inside the body
   j = logical((xx >= 0) .* (xx <= 1) .* (yy >= 0) .* (yy <= 1));
   xx = xx(j);
   yy = yy(j);
end

