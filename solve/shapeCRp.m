function val = shapeCRp(xp,c4e)
    v1 = (c4e(2,:) - c4e(1,:))';
    v2 = (c4e(3,:) - c4e(1,:))';
    xip = [v1,v2] \ (xp - c4e(1,:))';
    D = ones(3,3) - 2*eye(3);
    val0 = D*[1-sum(xip); xip];
    val = (val0+1)/4;
end