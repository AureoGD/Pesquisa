function [ Regions ] = partition_rest_regions( Rest_regions, G, W, S, H, E, F, tol, Nu, old_out_region, max )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    Regions = {};
    for i = 1:size(Rest_regions,1)
        A = Rest_regions{i,1};
        b = Rest_regions{i,2};
        [ xc , r, diagnostics] = chebychev_ball( A, b, G, W, S, H, F );  %Find epsilon and x0
        %plot(xc(1),xc(2),'*')
        
        if ((r <= 0.1))
            %figure(4)
            %hold on
             %Regions = [Regions; {Rest_regions{i,:}} max];
        elseif (isnan(r) || isnan(xc(1)) || (diagnostics.problem==1))
             %Regions = [Regions; {Rest_regions{i,:}} max];
             %plotregion(-A, -b)
        %elseif max>2    
         %    Regions = [Regions; {Rest_regions{i,:}} 6];
        else
            [ z0, diagnostics ] = optimal_z_mp_QP( G, W, S, H, F, xc, Nu);  % Find z0
            %z0 = fcnKKT(H, F, G, E, W, xc)
            
            G_tio = [];
            %if isempty(z0) == 0
            if diagnostics.problem == 0
                [G_tio W_tio S_tio] = verify_active_constraints(G, W, S, xc, z0, tol); %Find G_tio, W_tio e S_tio
            end
            %Define U(x)
            if isempty(G_tio) == 1
                Regions = [Regions; {Rest_regions{i,:}} max];
                x_lqr = xc;
                %plotregion(-A, -b)
                    
            elseif rank(G_tio) < size(G_tio,1)
                Regions = [Regions; {Rest_regions{i,:}} 20];
            else
                [ A, b ] = define_region( G, W, S, G_tio, W_tio, S_tio, H, tol); %Define polyhedron
                [A, b] = remove_redundant_constraints(A,b);
                CR = {A b};
                Regions = [Regions; CR 10];
                
                %plotregion(-A, -b)
                
                %Rest_CR = new_rest_regions(A,b,{Rest_regions{i,:}}, old_out_region);
                Rest_CR = find_rest_regions(A,b,{Rest_regions{i,:}});
              
                new_regions = partition_rest_regions(Rest_CR, G, W, S, H, E, F, tol, Nu, {Rest_regions{i,:}}, max+1);
                Regions = [Regions; new_regions];
            end
        end
    end
end

