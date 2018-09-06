function [ z0, lambda0 ] = fcnKKT(H, F, G, E, W, x)
%UNTITLED Summary of this function goes here
%   This is a Karush Kuhn Tucker function that calcules the optimun z of
%   system, respecting theirs restrictions.
%   Imput arguments:
%                   [H]
%                   [F]
%                   [G]
%                   [E]
%                   [W]
%                   [x]
%   Output argumets:
%                   [Zo]
%                   [Lambda0]
z0 = [];
lambda0 = [];
S = E+G*inv(H)*F';
nG = size(G);

for i=1:nG
    lambda(i,1)=-inv(G(i,:)*inv(H)*G(i,:)')*(W(i,:)+S(i,:)*x);
end

for i=1:nG
    z(:,i)= -inv(H)*G(i,:)'*lambda(i);
end

flag = 0;
for i=1:nG
    if lambda(i) >= 0
        aux = G*z(:,i)-W-S*x
        %Gz< W +Sx    Gz - W - Sx < 0
        for j=1:nG
            %test = aux(j)
            if aux(j) > 0
                flag = 1;
                break
            end
        end
        if flag == 0
            lambda0 = lambda(i);
            z0 = z(:,i);
        end
    end
end

end

