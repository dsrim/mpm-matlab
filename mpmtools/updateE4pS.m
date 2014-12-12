function e4pn = updateE4pS(c4n,n4e,e4n,x4p,e4p)
% update which element each material pt belongs to.
nrElems = size(n4e,1);
i4elem = 1:nrElems;
nrPts = size(x4p,1);
for j = 1:nrPts
  if ispine(x4p(j,:),c4n(n4e(e4p(j,:),:),:))
    e4pn(j,e4p(j,:)) = true;
  else
    c4e = [];
    e2search = unique(e4n(n4e(e4p(j,:),:),:));
    e2search = e2search(e2search~=0);
    for k = 1:length(e2search)
      c4e(:,:,k) = c4n(n4e(e2search(k),:),:);
    end
    elem = e2search(ispine(x4p(j,:),c4e));
    e4pn(j,elem) = true;
    if isempty(elem)
      display('empty row in e4p')
    end
  end
end
end