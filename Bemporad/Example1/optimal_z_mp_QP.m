function [ z0, diagnostics ] = optimal_z_mp_QP( G, W, S, H, F, x0, Nu)
%Calculate the optmizal solution to the mp-QP problem: 
%V(x) = min 0.5*z'*H*z , subject to Gz <= W +S*x(t)
    U = sdpvar(Nu,1,'full'); 

    z = U + inv(H)*F'*x0; 

    LMI = [];
    %LMI = [LMI, -G*z + W + S*x0 >= 0 ];
     LMI = [LMI, G*z <= W + S*x0 ];
    objetivo = 0.5*z'*H*z;
    options = sdpsettings;
    options.solver = 'sedumi';

    diagnostics = optimize(LMI,objetivo,options)
    z0 = double(z)

end

