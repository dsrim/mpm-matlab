function val = ispine(p,c4e)
% check if p belongs to an element of given coordinates.
e2chk = size(c4e,3);
val = logical(zeros(e2chk,1)); 
for j = 1:e2chk
  vecs = c4e(:,:,j) - repmat(p,4,1);
  nvecs = vecs./repmat(sqrt(sum(vecs.^2,2)),1,2);
  angles = nvecs*nvecs';
  if sum(abs(acos([diag(angles,1); angles(4,1)])),1) >= (2*pi - 1e-8)
    val(j) = true;
    break
  end
end
end