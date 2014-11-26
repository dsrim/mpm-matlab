function val = shapeRN0(xp,c4n,nodes,bdNormal)
  h = [norm(c4n(4,:) - c4n(1,:)), norm(c4n(2,:) - c4n(1,:))];
  xi = (xp - c4n(1,:))./h; 
  loop = nodes;
  bdpart = [];
  for j=1:4
    k = (nodes(j) == bdNormal);
    if ~isempty(bdpart) && sum(k) > 0
      bdpart = [bdpart; bdNormal(k(:,1),:)];
    end
  end
  val = [(1-xi(1))*(1 - xi(2)); xi(1)*(1-xi(2));
         xi(1)*xi(2);(1-xi(1))*xi(2)];
  for j = 1:length(bdpart)
    val(bdpart(j,1),bdpart(j,2)) = 0;
  end
       
end