function [ G_tio, W_tio, S_tio ] = verify_active_constraints(G, W, S, x0, z0, tol)
%Find which constraints are active at given points x0 and z0

    index = [];
    for i = 1:length(W)
        (G(i,:)*z0 - W(i,:) - S(i,:)*x0)
        if ((G(i,:)*z0 - W(i,:) - S(i,:)*x0 < tol) && (double(G(i,:)*z0 - W(i,:) - S(i,:)*x0)> -tol))
            index = [index , i];
        end
    end
    G_tio = [];
    S_tio = [];
    W_tio = [];
    if (length(index)>0)
        for i = 1:length(index)
            G_tio = [G_tio; G(index(i),:)];
            S_tio = [S_tio; S(index(i),:)];
            W_tio = [W_tio; W(index(i),:)];
        end
    end
end

