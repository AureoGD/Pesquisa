function [ CRest ] = find_rest_regions( A, b, out_region )
%UNTITLED2 Summary of this function goes here
  %Detailed explanation goes here
    A_Ri_rest = [];
    b_Ri_rest = [];
    %CRest = {};
    Nx = 0;
    for i=1:(size(A,1))
        if (A(i,:) == [1 0])
            Nx= Nx+1; 
        elseif (A(i,:) == [0 1])
            Nx= Nx+1; 
        elseif (A(i,:) == [-1 0])
            Nx= Nx+1; 
        elseif (A(i,:) == [0 -1])  
            Nx= Nx+1;     
        end
    end

    Nr=Nx
    for i = 1:(size(A,1)-Nr)
         A_Ri_rest = [];
         b_Ri_rest = [];
            if (size(A,1)-Nr)>1
                for j = 1:(i-1)
                    %CRest{i,1} = [CRest{i,1}; A(j,:)];
                    %CRest{i,2} = [CRest{i,2}; b(j,:)];
                    A_Ri_rest = [A_Ri_rest ; A(j,:)];
                    b_Ri_rest = [b_Ri_rest ; b(j,:)];
                end
            end
            i
            %A(i,:)
            %CRest{i,1} = [CRest{i,1}; -A(i,:)];
            %CRest{i,2} = [CRest{i,2}; -b(i,:)];
            A_Ri_rest = [A_Ri_rest ; -A(i,:)]
            b_Ri_rest = [b_Ri_rest ; -b(i,:)]
            
            new_A =  A_Ri_rest;
            new_b =  b_Ri_rest;
            
            if isempty(out_region) == 0
                %new_A = CRest{i,1};
                %new_A = new_A(1:(size(new_A,1)-Nx),:);
  
                
                %new_b = CRest{i,2};
                %new_b = new_b(1:(size(new_b,1)-Nx),:);

                new_A = [new_A; out_region{1,1}];
                new_b = [new_b; out_region{1,2}];
                %CRest{i,1} = [new_A; A((size(A,1)-Nx+1):size(A,1),:); out_region{1,1}];
                %CRest{i,2} = [new_b; b((size(b,1)-Nx+1):size(b,1),:); out_region{1,2}];

            end
            [CRest{i,1}, CRest{i,2}] = remove_redundant_constraints(new_A,new_b);
            
    end
    


end

