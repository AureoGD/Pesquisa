function [ xc , r, diagnostics] = chebychev_ball( A, b, G, W, S, H, F )
%Return the center of the largest possible ball that can be placed inside
%the region defined by Ax<=b.
    

    xc = sdpvar(2,1,'full');
    r = sdpvar(1);
    %U = sdpvar(2,1,'full'); 
    z = sdpvar(2,1,'full'); 
    %z = U + inv(H)*F'*xc; 
    LMI = [];
    for i = 1:size(A,1)
        LMI = [LMI; A(i,:)*xc + r*sqrt(sum(A(i,:).^2)) <= b(i)];
    end
    LMI = [LMI; G*z - S*xc <= W];
%     optimize(A*xc+r*sqrt(sum(A.^2,2)) <= b,-r)

    options=sdpsettings;
    options.solver='sedumi';
    diagnostics = optimize(LMI,-r,options);
    xc = double(xc);
    r = double(r);

end

