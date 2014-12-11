function allsides = bdNormals4s(bdNormals,s4n,s4e,e4s)
    cpnts = [-2,1,2,-1];
    allsides = [];
    for j = cpnts
        elems = [];
        ind = bdNormals((bdNormals(:,2) == j),1);
        nodes = [ind(1:(end-1)) ind(2:end)];
        for k = 1:size(nodes,1)
            elems = [elems; e4s(s4n(nodes(k,1),nodes(k,2)))];
        end
        sides = s4e(elems,:); 
        allsides = [allsides; sides(:), j*ones(length(sides(:)),1)];
    end
    

end
% 
%         for k = 1:size(nodes,1)
%            sides = s4n(nodes(k,1),nodes(k,2));
%            allsides = [allsides; sides(:),j*ones(length(sides(:)),1)];
%         end

