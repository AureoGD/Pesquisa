function [ xc , r] = chebychev_ball( A, b )
%Return the center of the largest possible ball that can be placed inside
%the region defined by Ax<=b.

    xc = sdpvar(2,1);
    r = sdpvar(1);
    optimize(A*xc+r*sqrt(sum(A.^2,2)) <= b,-r)
    xc = double(xc);
    r = double(r);

end

