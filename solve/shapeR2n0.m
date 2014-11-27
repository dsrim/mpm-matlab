function val = shapeR2n0(xp,c4nL,nodes,bdNormals)
  h = [norm(c4nL(4,:) - c4nL(1,:)), norm(c4nL(2,:) - c4nL(1,:))];
  xi = (xp - c4nL(1,:))./h; 
  val = [(1-xi(1))*(1 - xi(2));
         xi(1)*(1-xi(2));
         xi(1)*xi(2);
         (1-xi(1))*xi(2)];
  val = repmat(val,1,2);
  for j = 1:4
    ind = (nodes(j) == bdNormals(:,1));
    val(j,abs(bdNormals(ind,2))) = 0;
  end
end