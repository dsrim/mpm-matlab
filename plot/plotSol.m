function plotSol(c4n,n4e,x4p,t,v4p,vpexact,vI,vIexact)

nrElems = size(n4e,1);
figure(2);clf;
for elem = 1:nrElems
  plot(c4n(n4e(elem,[1:4 1]),1),c4n(n4e(elem,[1:4,1]),2));
  hold on;
end
plot(x4p(:,1),x4p(:,2),'r.')
axis equal
margin = .01;
axis([min(c4n(:,1))-margin,max(c4n(:,1))+margin,...
      min(c4n(:,2))-margin,max(c4n(:,2))+margin])
title({'Mesh';['t=' num2str(t,'%10.3f')]})
% quiver(x4p(:,1),x4p(:,2),v4p(:,1),v4p(:,2),'r');
% quiver(x4p(:,1),x4p(:,2),vpexact(:,1),vpexact(:,2),'g');
quiver(c4n(:,1),c4n(:,2),vI(:,1),vI(:,2),1.5,'r');
quiver(c4n(:,1),c4n(:,2),vIexact(:,1),vIexact(:,2),1.5,'g');

end