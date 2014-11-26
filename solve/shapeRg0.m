function val = shapeRg0(xp,c4n,nodes,bdNodes)
  h = [norm(c4n(2,:) - c4n(1,:)), norm(c4n(4,:) - c4n(1,:))];
  xi = (xp - c4n(1,:))./h; 
  loop = nodes;
  bdpart = [];
  for j=1:4
    if ~isempty(bdpart) && sum(bdpart(j) == bdNodes) > 0
      bdpart = [bdpart; j];
    end
  end
  val = [ -(1-xi(2))/h(1), -(1-xi(1))/h(2);
          (1-xi(2))/h(1),   -xi(1)/h(2);
          xi(2)/h(1), xi(1)/h(2); 
          -xi(2)/h(1), (1-xi(1))/h(2)];
  val(bdpart,:) = zeros(length(bdpart),2);
       
end