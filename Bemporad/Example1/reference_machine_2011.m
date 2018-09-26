clear all
close all
clc

%SetPoint Control
%System from Example 1 Bemporad 2002
Ap = [0.7326 -0.0861;
     0.1722 0.9909];
Bp = [0.0609; 0.0064];
Cp = [0 1.4142];

A = [Ap zeros(1,2)';
     Cp*Ap eye(1)]
B = [Bp;
    Cp*Bp]
C = [zeros(1,2) eye(1)];
 
 %%
Ny = 3;
Nu = 3;
Nc = 3;
Nstate = 5;

Q = eye(1);
R = 0.01;

Umax = 60;
Umin = -60;
delta_Umax = 30;
delta_Umin = -30;
ref_max = 4;
ref_min = -4;

x0 = [1 1];
x_inicial = [x0 Cp*x0' 0 0.7]';

phi = [];
for i = 1:Ny
    phi = [phi; C*A^i];
end

%%

theta = [];

for i = 1:Ny
    for j = 1:Nu
       if j-i ==0
            theta(i,j) = C*B;
       elseif j < i
           theta(i,j) = C*A^(i-j)*B;
       else
           theta(i,j) = 0;
       end
        
    end
end

%%
Rlinha = [];
for i = 1:Nu
    Rlinha = blkdiag(Rlinha,R);
end

Qlinha = [];
for i = 1:Ny
    Qlinha = blkdiag(Qlinha,Q);
end
%%
Qt = ones(1,Nu);
%Qp = ones(

% lambda1 = [];
% lambda2 = [];
% for i = 1:(Ny-1+2)
%     lambda1 = [lambda1 (-i+2)*eye(1)];
%     lambda2 = [lambda2 (i-1)*eye(1)];
% end
% lambda = [lambda1; lambda2]
% %thetaP = theta((Ny-1+2):size(theta,1),:)
% hue = -ones(1,Ny)*theta;

H = 2*Rlinha + 2*theta'*Qlinha*theta;

F = [2*phi'*Qlinha*theta;
    zeros(1,Nu);
    -2*Qt*theta];
    
%Restrição deltaU
G = [1;
     -1];
W = [delta_Umax;
      -delta_Umin];
E = zeros(2*Nu,5); 
 for i = 1:(Nu-1)
    G = blkdiag(G,[1; -1]);
    W = [W; delta_Umax; -delta_Umin];
 end
 
%Restricao U
G = [G ; tril(ones(Nu));- tril(ones(Nu))];
for i = 1:Nu
    E = [E; 0 0 0 -1 0];
    W = [W; Umax];
end
for i = 1:Nu 
    E = [E; 0 0 0 1 0];
    W = [W; -Umin];
end

%Restrição u(t-1)
G = [G; zeros(2*Nu,Nu)];
out_Aold = [];
out_Bold = [];
for i = 1:Nu
    W = [W; -Umin ; Umax];
    E = [E; 0 0 0 1 0; 0 0 0 -1 0];
    out_Aold = [out_Aold; 0 0 0 -1 0; 0 0 0 1 0];
    out_Bold = [out_Bold; -Umin ; Umax];
end

%Restrição r(t)
G = [G; zeros(2*Nu,Nu)];
out_Ar = [];
out_Br = [];
for i = 1:Nu
    W = [W; -ref_min ; ref_max];
    E = [E; 0 0 0 0 1; 0 0 0 0 -1];
    out_Ar = [out_Ar; 0 0 0 0 -1; 0 0 0 0 1];
    out_Br = [out_Br; -ref_min ; ref_max];
end

S = E + G*inv(H)*F'; 

%%
z0 = optimal_z_mp_QP( G, W, S, H, F, x_inicial, Nu);
%%
tol = 1e-5;
[G_tio, W_tio, S_tio] = verify_active_constraints(G, W, S, x_inicial, z0, tol);
%%
[A_CR0 b_CR0] = define_region(G, W, S, G_tio, W_tio, S_tio, H, tol);
%%
[A_CR0 b_CR0] = remove_redundant_constraints(A_CR0, b_CR0, Nu, Nstate)
%%
[Kx, Kc] = define_control(G, W, S, G_tio, W_tio, S_tio, H, F);
Regions{1,1} = A_CR0;
Regions{1,2} = b_CR0;
Regions{1,3} = Kx;
Regions{1,4} = Kc;
%%
out_X{1,1} = [out_Aold; out_Ar];
out_X{1,2} = [out_Bold; out_Br];

%%
CR0_rest = find_rest_regions(A_CR0, b_CR0, Nu, Nstate, out_X);

%%
new_Regions = partition_rest_regions(CR0_rest, G, W, S, H, F, tol, Nu, Nstate, 1);
Regions = [Regions; new_Regions];
