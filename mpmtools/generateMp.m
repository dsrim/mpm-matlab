function [xp,e4p] = generateMp(c4n,n4e,isbody,nmpe)
mmpe = ceil(sqrt(nmpe));
e4p = logical([]);
xp = [];
nrElems = size(n4e,1);
for elem = 1:nrElems
  mp4e = [];
  xlim4e = c4n(n4e(elem,[1,2]),1);
  ylim4e = c4n(n4e(elem,[2,3]),2);
  xx = linspace(xlim4e(1),xlim4e(2),mmpe+2);
  yy = linspace(ylim4e(1),ylim4e(2),mmpe+2);
  [X,Y] =  meshgrid(xx(2:(end-1)),yy(2:(end-1)));
  [mp4e(:,1),mp4e(:,2)] = isbody(X(:),Y(:));
  xp = [xp; mp4e];
  k = size(xp,1);
  npts = size(mp4e,1);
  e4p((k-npts+1):k,elem) = ones(npts,1);
end

end