function plotSol(c4n,n4e,x4p,t,v4p,wI)
% 
nrElems = size(n4e,1);
nrPts = size(x4p,1);
f = figure('visible','off'); 
clf;
h = max(norm(c4n(3,:) - c4n(1,:),2));
scale = 60*h;
for elem = 1:nrElems
  plot(c4n(n4e(elem,[1:4 1]),1),c4n(n4e(elem,[1:4,1]),2));
  hold on;
end
axis equal
margin = .01;
axis([min(c4n(:,1))-margin,max(c4n(:,1))+margin,min(c4n(:,2))-margin,max(c4n(:,2))+margin])
title({'MPM Solution u(X,t)'; ['t=' num2str(t,'%10.3f') ',' num2str(nrElems) ' elts,' num2str(nrPts) ' mtrl pts']})

quiver(c4n(:,1),c4n(:,2),wI(:,1),wI(:,2),h,'r');
hold off;
legend('mesh','u_I','Location','NorthEastOutside');
drawnow
% matlab2tikz(['../output/mpm-sol-' num2str(nrElems) '-' num2str(nrPts/nrElems) ...
%     '-' num2str(t,'%1.4f') '.tikz'])

f = figure('visible','off'); 
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
quiver(x4p(:,1),x4p(:,2),v4p(:,1),v4p(:,2),'r');

title({'MPM Solution v(x_p,t)';['t=' num2str(t,'%10.3f') ',' num2str(nrElems) ' elts,'...
   num2str(nrPts) ' mtrl pts']})
legend('mesh','mtrl pts','v_I','Location','NorthEastOutside');
drawnow
%matlab2tikz(['../output/mpm-dsoldt-' num2str(nrElems) '-' num2str(nrPts/nrElems) ...
    %'-' '-' num2str(t,'%1.4f') '.tikz'])
saveas(f,['../output/mpm-dsoldt-' num2str(nrElems) '-' num2str(nrPts/nrElems) '-' num2str(t,'%1.4f') '.eps'], 'eps')

f = figure('visible','off'); clf;
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
% matlab2tikz(['../output/mpm-grid-' num2str(nrElems) '-' num2str(nrPts/nrElems) ...
%     '-' '-' num2str(t,'%1.4f') '.tikz'])
end
