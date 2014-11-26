
c4n = [0,0;1,0;1,1;0,1;];
target = [4,0;5,4;3,6;2,4];
xp = mean(target,1);

figure(1);
plot(c4n([1:4,1],1),c4n([1:4,1],2),'-bo',target([1:4,1],1),target([1:4,1],2),'-rx',xp(1),xp(2),'bx');
xi = linspace(0,1,5);
[Xi1,Xi2] = meshgrid(xi,xi);
figure(2);
b = (xp - target(1,:))';
A = [(target(2,:) - target(1,:))', ...
     (target(4,:) - target(1,:))', ...
     (target(3,:) - target(1,:))'];
xi = A(:,[1 2]) \ b;
xi1 = xi(1);

b = b - A(:,2)*xi1;
B = A(:,2)+A(:,3)*xi1;

