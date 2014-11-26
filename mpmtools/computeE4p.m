function ne4p = computeE4p(c4n,n4e,e4n,xp,e4p)

nrElems = size(n4e,1);
i4elem = 1:nrElems;
nrPts = size(xp,1);
for j = 1:nrPts
  if ispine(xp(j,:),c4n(n4e(e4p(j,:),:),:))
    ne4p(j,e4p(j,:)) = true;
  else
    c4e = [];
    e2search = i4elem(sum(e4n(n4e(e4p(j,:),:),:),1) > 0);
    for k = 1:length(e2search)
      c4e(:,:,k) = c4n(n4e(e2search(k),:),:);
    end
    elem = e2search(ispine(xp(j,:),c4e));
    ne4p(j,elem) = true;
    if isempty(elem)
      display('empty row in e4p')
    end
  end
end


end