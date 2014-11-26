function val = shape(xp,c4e)
  xi = (xp - c4e(1,:)); % needs more work!
  val = [(1-xi(1))*(1 - xi(2));
         xi(1)*(1-xi(2));
         xi(1)*xi(2);
         (1-xi(1))*xi(2)];
end