function plotSol(c4n,n4e,x4p,t,v4p,vpexact,wI,wIexact)

nrElems = size(n4e,1);
nrPts = size(x4p,1);
figure(1);clf;
h = max(norm(c4n(3,:) - c4n(1,:),2));
scale = 60*h;
for elem = 1:nrElems
  plot(c4n(n4e(elem,[1:4 1]),1),c4n(n4e(elem,[1:4,1]),2));
  hold on;
end
axis equal
margin = .01;
axis([min(c4n(:,1))-margin,max(c4n(:,1))+margin,...
      min(c4n(:,2))-margin,max(c4n(:,2))+margin])
title({'MPM Solution u(X,t)';['t=' num2str(t,'%10.3f') ',' num2str(nrElems) ' elts,'...
    num2str(nrPts) ' mtrl pts']})

quiver(c4n(:,1),c4n(:,2),scale*wIexact(:,1),scale*wIexact(:,2),0,'--g');
quiver(c4n(:,1),c4n(:,2),scale*wI(:,1),scale*wI(:,2),0,'r');
hold off;
legend('mesh','u_I','uexact','Location','NorthEastOutside');
drawnow
matlab2tikz(['output/mpm-sol-' num2str(nrElems) '-' num2str(nrPts/nrElems) ...
    '-' num2str(t,'%1.4f') '.tikz'])
% quiver(c4n(:,1),c4n(:,2),10*(wI(:,1)-wIexact(:,1)),10*(wI(:,2)-wIexact(:,2)),'r');

figure(2); 
clf;
for elem = 1:nrElems
  plot(c4n(n4e(elem,[1:4 1]),1),c4n(n4e(elem,[1:4,1]),2));
  hold on;
end
plot(x4p(:,1),x4p(:,2),'r.');

axis equal
margin = .01;
axis([min(c4n(:,1))-margin,max(c4n(:,1))+margin,...
      min(c4n(:,2))-margin,max(c4n(:,2))+margin])
title({'MPM Solution';['t=' num2str(t,'%10.3f')]})
quiver(x4p(:,1),x4p(:,2),vpexact(:,1),vpexact(:,2),'--g');
quiver(x4p(:,1),x4p(:,2),v4p(:,1),v4p(:,2),'r');

title({'MPM Solution v(x_p,t)';['t=' num2str(t,'%10.3f') ',' num2str(nrElems) ' elts,'...
    num2str(nrPts) ' mtrl pts']})
legend('mesh','mtrl pts','v_I','vexact','Location','NorthEastOutside');
drawnow
matlab2tikz(['output/mpm-dsoldt-' num2str(nrElems) '-' num2str(nrPts/nrElems) ...
    '-' '-' num2str(t,'%1.4f') '.tikz'])

figure(3); clf;
h = max(norm(c4n(3,:) - c4n(1,:),2));
for elem = 1:nrElems
  plot(c4n(n4e(elem,[1:4 1]),1)+wI(n4e(elem,[1:4,1]),1),...
      c4n(n4e(elem,[1:4,1]),2)+wI(n4e(elem,[1:4,1]),2));
  hold on;
end
axis equal
margin = .01;
axis([min(c4n(:,1))-margin,max(c4n(:,1))+margin,...
      min(c4n(:,2))-margin,max(c4n(:,2))+margin])
title({'MPM Solution u(x,t)';['t=' num2str(t,'%10.3f') ',' num2str(nrElems) ' elts,'...
    num2str(nrPts) ' mtrl pts']})
drawnow
matlab2tikz(['output/mpm-grid-' num2str(nrElems) '-' num2str(nrPts/nrElems) ...
    '-' '-' num2str(t,'%1.4f') '.tikz'])
end