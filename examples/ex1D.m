function ex1D

E = 1; 
L = 1;
A = 0.01;
v0 = @(x) A*sin(pi*x);
u0 = @(x) 0*x;
rho0 = @(x) 0*x+1;
b0 = @(x) 0*x;
sigma0 = @(x) 0*x ;
ppe = 2;
nrElems = 8;
nrNodes = nrElems+1;
nrSteps = 20;
T = 10;
dt = T/nrSteps;
c4n = linspace(0,L,nrNodes)';
h = L/(nrNodes-1);
n4e = [1:nrElems; 2:(nrElems+1)]';
e4n = sparse(n4e(:,1),n4e(:,2),1:nrElems);
bdNodes = [1 nrNodes];
inNodes = setdiff(1:nrNodes,bdNodes);
nrPts = nrElems*ppe;

%% Initialization
xp = zeros(nrPts,nrSteps);
e4xp = zeros(nrPts,nrSteps);
volp = zeros(nrPts,nrSteps);
rhop = zeros(nrPts,nrSteps);
xp = zeros(nrPts,nrSteps);
vp = zeros(nrPts,nrSteps);
bp = zeros(nrPts,nrSteps);
sigmap = zeros(nrPts,nrSteps);
Fp = zeros(nrPts,nrSteps);
Jp = zeros(nrPts,nrSteps);
FI = zeros(nrNodes,nrSteps);
vI = zeros(nrNodes,nrSteps);
uI = zeros(nrNodes,nrSteps);
xpt = zeros(ppe,nrElems); volpt = zeros(ppe,nrElems);
icon = diag(ones(nrElems-1,1),-1) + diag(ones(nrElems-1,1),1);
for elem = 1:nrElems
   a = c4n(n4e(elem,1)); b = c4n(n4e(elem,2));
   pts = linspace(a,b,ppe+2); 
   xpt(:,elem) = pts(2:(end-1))';
   xp4et(elem,:) = elem;
   volpt(:,elem) = abs(a-b)/2;
end
xp(:,1) = xpt(:); 
xp0 = xp(:,1);
xp4e(:,1) = xp4et(:);
volp(:,1) = volpt(:);
bp(:,1) = b0(xp(:,1));
vp(:,1) = v0(xp(:,1));
up(:,1) = u0(xp(:,1));
rhop(:,1) = rho0(xp(:,1));
rhop0 = rhop(:,1);
sigmap(:,1) = sigma0(xp(:,1));
mp(:,1) = volp(:,1).*rhop(:,1);
Fp(:,1) = ones(nrPts,1);
Jp(:,1) = abs(Fp(:,1));
e4xp0 = gete4xp(xp0,c4n,e4n);

%% MPM Algorithm


for n = 1:(nrSteps-1)
    % Map to grid
    e4xp = gete4xp(xp(:,n),c4n,e4n);
    [mI,vI(:,n),fiI,fxI,fI,FI(:,n+1)] = p2I(xp(:,n),e4xp,c4n,n4e,mp(:,1),volp(:,n), ...
                            vp(:,n),bp(:,n),sigmap(:,n),Fp(:,n),...
                            nrPts,nrNodes,inNodes,h);    
    % Solve / Timestep
    vI(1,n) = 0;
    vI(end,n) = 0;
    vI(inNodes,n+1) = vI(inNodes,n) + dt*fI(inNodes)./mI(inNodes);
    uI(inNodes,n+1) = uI(inNodes,n) + dt*vI(inNodes,n+1);

    % Update material points
    [xp(:,n+1),vp(:,n+1),Fp(:,n+1),Jp(:,n+1),...
        rhop(:,n+1),volp(:,n+1),sigmap(:,n+1)] = updatePts(...
                        xp0,xp,vp(:,n),mp,mI,volp(:,n),e4xp,e4xp0,c4n,n4e,rhop0,vI(:,[n n+1]),Fp(:,n),...
                        FI(:,n+1),nrPts,nrNodes,h,dt);
    xx = xp(:,n+1);
%     plot(xx,vp(:,n+1),xx,rhop(:,n+1),c4n,vI(:,n+1),c4n,FI(:,n+1))
    figure(1); plot(c4n,uI(:,n+1)); axis([0 1 0 1])
end





end

function [mI,vI,fiI,fxI,fI,FI] = p2I(xp,e4xp,c4n,n4e,mp,volp,vp,...
                                            bp,sigmap,Fp,nrPts,nrNodes,inNodes,h)
mI = zeros(nrNodes,1);
vI = zeros(nrNodes,1);
fiI = zeros(nrNodes,1);
fxI = zeros(nrNodes,1);
FI = zeros(nrNodes,1);
for p = 1:nrPts
    Ns = n4e(e4xp(p),:);
    Nxp = shape(xp(p),Ns,c4n,h);
    dNxpdx = shapeg(xp(p),Ns,c4n,h);
    NxpAll(:,:,p) = Nxp;
    onecol = ones(length(Ns),1);
    mIloc = mp(p).*shapef(xp(p),Ns,c4n,h);
    mI = mI + sparse(Ns, onecol, mIloc,nrNodes,1);
    vIloc = vp(p).*mp(p).*Nxp;
    vI = vI + sparse(Ns,onecol, vIloc,nrNodes,1);
    fiIloc = sigmap(p).*volp(p).*dNxpdx;
    fiI = fiI - sparse(Ns,onecol,fiIloc,nrNodes,1);
    fxIloc = mp(p).*bp(p).*shapef(xp(p),Ns,c4n,h);
    fxI = fxI + sparse(Ns,onecol,fxIloc,nrNodes,1);
%     FIloc = volp(p).*Fp(p).*shapef(xp(p),Ns,c4n,h);
%     FI = FI + sparse(Ns,onecol,FIloc,nrNodes,1);
end
vI = vI./mI;
% FI = FI./mI;
fI = fiI + fxI;
 end

function [xpn,vpn,Fpn,Jpn,rhopn,volpn,sigmapn] = updatePts(xp0,xp,vp,mp,mI,volp,e4xp,e4xp0,c4n,n4e,...
                            rhop0,vI,Fp,FI,nrPts,nrNodes,h,dt)

% needs work
E = 1;
Fpfactor = zeros(nrPts,1);
xpn = zeros(nrPts,1);
vpn = zeros(nrPts,1);
sigmafactor = zeros(nrPts,1);
for p = 1:nrPts
    Ns = n4e(e4xp(p),:);
    NIxp = shape(xp(p),Ns,c4n,h);
    dNIxp = shapeg(xp(p),Ns,c4n,h);
    Fpfactor(p) = 1 + dt*sum(vI(Ns,2).*shapegf(xp(p),Ns,c4n,h));
    xpn(p) = xp(p) + dt*sum(vI(Ns,2).*NIxp);
    vpn(p) = vp(p) + sum((vI(Ns,2)-vI(Ns,1)).*NIxp);
    FIloc = volp(p).*Fp(p).*shapef(xp(p),Ns,c4n,h);
    FI = FI + sparse(Ns,ones(length(Ns),1),FIloc,nrNodes,1);
end

for p = 1:nrPts
    Ns0 = n4e(e4xp0(p),:);
    sigmafactor(p) = sum(FI(Ns0).*shapef(xp0(p),Ns0,c4n,h));
end
Fpn = Fpfactor.*Fp;
Jpn = abs(Fpn);
rhopn = rhop0./Jpn;
volpn = mp./Jpn;
sigmapn = E.*sigmafactor;

end

function y = shape(xp, Ns, c4n, h)
    y = 1- abs(xp - c4n(Ns))/h;
    if sum(Ns == 1)
       y(Ns ==1 ) = 0;
    elseif sum(Ns == length(c4n))
        y(Ns == length(c4n)) = 0;
    end
end

function y = shapef(xp, Ns, c4n, h)
    y = 1- abs(xp - c4n(Ns))/h;

end

function y = shapeg(xp, Ns, c4n, h)
    y = -sign(xp - c4n(Ns))/h;
    if sum(Ns == 1)
       y(Ns ==1 ) = 0;
    elseif sum(Ns == length(c4n))
        y(Ns == length(c4n)) = 0;
    end
end

function y = shapegf(xp, Ns, c4n, h)
    y = -sign(xp - c4n(Ns))/h;

end

function val = gete4xp(xp,c4n,e4n)
    nrNodes = length(c4n);
    nrPts = length(xp);
    nodes1 = sum(xp*ones(1,nrNodes) - ones(nrPts,1)*c4n'>0,2);
    val = diag(e4n(nodes1,nodes1+1));
end
