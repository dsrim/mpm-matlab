function val = shapeR(xp,c4n)
  h = [norm(c4n(4,:) - c4n(1,:)), norm(c4n(2,:) - c4n(1,:))];
  xi = (xp - c4n(1,:))./h; 
  val = [(1-xi(1))*(1 - xi(2));
         xi(1)*(1-xi(2));
         xi(1)*xi(2);
         (1-xi(1))*xi(2)];
end