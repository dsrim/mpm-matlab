function val = ispine(p,c4e)
% check if p belongs to an element of given coordinates.
e2chk = size(c4e,3);
val = logical(zeros(e2chk,1)); 
nCols = size(c4e,1);
for j = 1:e2chk
  vecs = c4e(:,:,j) - repmat(p,nCols,1);
  nvecs = vecs./repmat(sqrt(sum(vecs.^2,2)),1,2);
  angles = nvecs*nvecs';
  m = size(angles,1);
  if sum(abs(acos([diag(angles,1); angles(m,1)])),1) >= (2*pi - 1e-8)
    val(j) = true;
    break
  end
end
end