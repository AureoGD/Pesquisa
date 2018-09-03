function [ Regions ] = partition_rest_regions( Rest_regions, G, W, S, H, F, tol, Nu, old_out_region, max )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    Regions = {};
    for i = 1:size(Rest_regions,1)
        A = Rest_regions{i,1};
        b = Rest_regions{i,2};
        [ xc , r, diagnostics] = chebychev_ball( A, b, G, W, S, H, F );  %Find epsilon and x0
        plot(xc(1),xc(2),'*')
        
        if ((r <= 0.1))
            figure(4)
            hold on
             Regions = [Regions; {Rest_regions{i,:}} 3];
        elseif (isnan(r) || isnan(xc(1)) || (diagnostics.problem==1))
            Regions = [Regions; {Rest_regions{i,:}} 10];
             %plotregion(-A, -b)
        elseif max>2    
             Regions = [Regions; {Rest_regions{i,:}} 6];
        else
            [ z0, diagnostics ] = optimal_z_mp_QP( G, W, S, H, F, xc, Nu);  % Find z0
            
            G_tio = [];
            if diagnostics.problem == 0
                [G_tio W_tio S_tio] = verify_active_constraints(G, W, S, xc, z0, tol); %Find G_tio, W_tio e S_tio
                6
            end
            %Define U(x)
            if isempty(G_tio) == 1
                Regions = [Regions; {Rest_regions{i,:}} 0];
                7
                x_lqr = xc;
                %plotregion(-A, -b)
                    
            elseif rank(G_tio) < size(G_tio,1)
                12
                Regions = [Regions; {Rest_regions{i,:}} 2];
            else
                
                8
                [ A, b ] = define_region( G, W, S, G_tio, W_tio, S_tio, H, tol); %Define polyhedron
                CR = {A b};
                Regions = [Regions; CR 1];
                
                %plotregion(-A, -b)
                
                Rest_CR = new_rest_regions(A,b,{Rest_regions{i,:}}, old_out_region);
                Nr = size(Rest_regions{i,1},1);
                Rest_CR{:,:};
                new_regions = partition_rest_regions(Rest_CR, G, W, S, H, F, tol, Nu, {Rest_regions{i,:}}, max +1);
                Regions = [Regions; new_regions];
            end
        end
    end
end

