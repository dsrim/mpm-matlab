function val = computeL2error(uI,uexact,c4n,n4e,T)
val = 0;
nrElems = size(n4e,1);
for elem = 1:nrElems
    a = c4n(n4e(elem,1),1);
    b = c4n(n4e(elem,2),1);
    c = c4n(n4e(elem,2),2);
    d = c4n(n4e(elem,3),2);
    c4nL = c4n(n4e(elem,:),:);
    uIL= uI(n4e(elem,:),:);
    val = val + ...
      quad2d(@(x,y) squareerror(x,y,uIL,c4nL,uexact,T),a,b,c,d,'AbsTol',1e-8);
end
val = sqrt(val);
end

function val = squareerror(x,y,uIL,c4nL,uexact,T)
for j = 1:size(x,1)
  for k = 1:size(x,2)
    val(j,k) = sum(abs(shapeR([x(j,k),y(j,k)],c4nL)'*uIL - uexact([x(j,k),y(j,k)],T)).^2);
  end
end

end