function e4n = computeE4nQ(n4e)
nrElems = size(n4e,1);
j = repmat(1:nrElems,1,4);
e4n = logical(sparse(n4e(:),j,ones(1,4*nrElems)));
end