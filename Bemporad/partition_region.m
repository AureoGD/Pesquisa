function [ Regions ] = partition_region( CRi_rest, G, W, S, H, F, Nu, Nx, tol, max, out_region )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    Regions = {};
    Nr = size(CRi_rest,1)
    for i = 1:Nr
        [xc , r] = chebychev_ball( CRi_rest{i,1}, CRi_rest{i,2} );
        xc;
        r;
        if (r > 0 && max < 3)
            [ z0 ] = optimal_z_mp_QP( G, W, S, H, F, xc, Nu);
            [ G_tio, W_tio, S_tio ] = verify_active_constraints(G, W, S, xc, z0, tol);
            [ A, b ] = define_region( G, W, S, G_tio, W_tio, S_tio, H, tol );
            CRest = rest_region( A, b, Nx, out_region)
            new_regions = partition_region(CRest, G, W, S, H, F, Nu, Nx, tol, max, {CRi_rest{i,:}});
            Regions = [Regions; new_regions];
        else
            Regions = [Regions; {CRi_rest{i,:}}];
        end
        2
    end
    3
end

