function plotMesh(c4n,n4e,x4p,t)

nrElems = size(n4e,1);
figure(1);clf;
for elem = 1:nrElems
  plot(c4n(n4e(elem,[1:4 1]),1),c4n(n4e(elem,[1:4,1]),2));
  hold on;
end
plot(x4p(:,1),x4p(:,2),'r.')
axis equal
margin = 0.01;
axis([min(c4n(:,1))-margin,max(c4n(:,1))+margin,...
      min(c4n(:,2))-margin,max(c4n(:,2))+margin])
title({'Mesh';['t=' num2str(t)]})

end