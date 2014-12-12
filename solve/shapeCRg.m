function grads = shapeCRg(c4e)
        grads = [c4e';1 1 1]\[-2 0; 0 -2; 0 0]; % gradients for CR basis
%         grads = 1/4*(ones(3,3) + eye(3))*grads;
        grads = grads/2;
end