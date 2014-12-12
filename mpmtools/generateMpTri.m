function [x4p,e4p] = generateMpTri(mid4s,s4e,isbody)
e4p = logical([]);
x4p = [];
nrElems = size(s4e,1);
for elem = 1:nrElems
  mp4e = [];
  mid4e = mid4s(s4e(elem,:),:);
  mp4e = (mid4e([1 2 3],:) + mid4e([2 3 1],:))/2;
%   mp4e = [mp4e; (2*mid4e + mid4e([2 3 1],:))/3];
  [mp4e(:,1),mp4e(:,2)] = isbody(mp4e(:,1),mp4e(:,2));
  x4p = [x4p; mp4e];
  k = size(x4p,1);
  npts = size(mp4e,1);
  e4p((k-npts+1):k,elem) = ones(npts,1);
end

end