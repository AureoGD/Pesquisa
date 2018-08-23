function [ CRest ] = partition_region( A, b, Nx )
%Partition CRest region from the CR0 matrices. Nx is the number of
%constraints that define the space.
    
    Nr = size(A,1) - Nx;
    %A_CRest{1:Nr,1} = A(1:(size(A,1)-Nr),:);
    
    %b_CRest =0;
    for i = 2:Nr
        %A_CRest{i} = A(1,:);
        CRest{i,1} = A(1,:);
        CRest{i,2} = b(1,:);
    end
    
    
    for i = 1:Nr    
        for j = 2:(i-1)
            %A_CRest{i} = [A_CRest{i}; A(j,:)];
            CRest{i,1} = [CRest{i,1}; A(j,:)];
            CRest{i,2} = [CRest{i,2}; b(j,:)];
        end
        %A_CRest{i} = [A_CRest{i}; -A(i,:); A((size(A,1)-Nx+1:size(A,1)),:)];
        CRest{i,1} = [CRest{i,1}; -A(i,:); A((size(A,1)-Nx+1:size(A,1)),:)];
        CRest{i,2} = [CRest{i,2}; -b(i,:); b((size(A,1)-Nx+1:size(A,1)),:)];
    end

end

