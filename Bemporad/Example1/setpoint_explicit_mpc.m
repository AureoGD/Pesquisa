clear all
close all
clc

%SetPoint Control
%System from Example 1 Bemporad 2002
A = [0.7326 -0.0861;
     0.1722 0.9909];
B = [0.0609; 0.0064];
C = [0 1.4142];

Ny = 5;
Nu = 5;
Nc = 5;
Nstate = 4;

Q = eye(Ny);
R = 0.01;

Umax = 2;
Umin = -2;
delta_Umax = 0.5;
delta_Umin = -0.5;
ref_max = 4;
ref_min = -4;

x0 = [1 1];
x_inicial = [x0 0 0.7]';
%P dado pela solução de Lyapunov Eq 4
%P = zeros(2);
%P = dlyap(A,Q)
%%
Qlinha = Q;%blkdiag(Q,1);
for i = 1:(Ny-1)
    C = blkdiag(C,[0 1.4142]);
end
%C = blkdiag(C,[0 1.4142])
% for i = 1:(Nu-1)
%     Qlinha = blkdiag(Qlinha,Q);
% end
% Qlinha = blkdiag(Qlinha,P);
Rlinha = R*eye(Nu);

%%
Sx = eye(length(A));
%for i = 1:Ny
for i = 1:(Ny-1)
    Sx = [Sx; A^i];
end

%Sd = zeros((Ny+1)*length(B),Nu);
Sd = zeros((Ny)*length(B),Nu);
for j = 1:Nu
    %for i = 1:(Ny+1)
    for i = 1:(Ny)
        if (i-j-1)>0
            Sd((length(B)*(i-1)+1):(length(B)*(i)),j) =  [A^(i-j-1)*B];
        elseif (i-j-1) == 0
            Sd((length(B)*(i-1)+1):(length(B)*(i)),j) = [B];
        end
    end
end

Su = sum(Sd,2);
%%
H = Sd'*C'*Qlinha*C*Sd + Rlinha;
%H = C*Sd'*Qlinha*Sd*C' + Rlinha;

%%
F = [Sx'*C'*Qlinha*C*Sd]
F = [F; Su'*C'*Qlinha*C*Sd]
F = [F; -ones(1,Ny)*Qlinha*C*Sd]
% F = [C*Sx'*Qlinha*Sd*C']
% F = [F; C*Su'*Qlinha*Sd*C']
% F = [F; -Qlinha*Sd*C']

  %%   
% G = [1 0;
%      -1 0;
%      0 1;
%      0 -1];
%  W = [delta_Umax;
%       -delta_Umin;
%       delta_Umax;
%       -delta_Umin];
% E = zeros(4,4);

%Restricao deltaU
G = [1;
     -1];
W = [delta_Umax;
      -delta_Umin];
E = zeros(2*Ny,4); 
 for i = 1:(Ny-1)
    G = blkdiag(G,[1; -1]);
    W = [W; delta_Umax; -delta_Umin];
 end

 %%
%Restricao U
G = [G ; tril(ones(Nu));- tril(ones(Nu))];

for i = 1:Nu
    E = [E; 0 0 -1 0];
    W = [W; Umax];
end
for i = 1:Nu 
    E = [E; 0 0 1 0];
    W = [W; -Umin];
end

%Restrição u(t-1)
G = [G; zeros(2*Nu,Nu)];
for i = 1:Nu
    W = [W; -Umin ; Umax];
    E = [E; 0 0 1 0; 0 0 -1 0];
end

%Restrição r(t)
G = [G; zeros(2*Nu,Nu)];
for i = 1:Nu
    W = [W; -ref_min ; ref_max];
    E = [E; 0 0 0 1; 0 0 0 -1];
end


%%
S = E + G*inv(H)*F'; 

%%
z0 = optimal_z_mp_QP( G, W, S, H, F, x_inicial, Nu)
%%
tol = 1e-8;
[G_tio, W_tio, S_tio] = verify_active_constraints(G, W, S, x_inicial, z0, tol);
%%
[A_CR0 b_CR0] = define_region(G, W, S, G_tio, W_tio, S_tio, H, tol);
%%
[Kx, Kc] = define_control(G, W, S, G_tio, W_tio, S_tio, H, F);
Regions{1,1} = A_CR0;
Regions{1,2} = b_CR0;
Regions{1,3} = Kx;
Regions{1,4} = Kc;
%out_X{1,1} = [1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0 ; 0 0 -1];
%out_X{1,1} = [1 0 ; 0 1 ; -1 0 ; 0 -1];
%out_X{1,2} = [Xmax' ; -Xmin'];
out_X = {};
%%
CR0_rest = find_rest_regions(A_CR0, b_CR0, Nu, Nstate, out_X);

%%
new_Regions = partition_rest_regions(CR0_rest, G, W, S, H, F, tol, Nu, Nstate, 1);
Regions = [Regions; new_Regions];

%%
figure(5)
hold on
for i=1:size(Regions,1)
    HUEA = Regions{i,1};
    HUEA = [HUEA(:,1:2); 1 0; -1 0; 0 1; 0 -1]
    HUEB = Regions{i,2};
    HUEB = [HUEB; 5 ; 5; 5; 5]
    plotregion(-HUEA,-HUEB)
    xlim([-5 5])
    ylim([-5 5])
end
%plotregion(-A_CR0(:,1:2),-b_CR0)
%xlim([-2 2])
%ylim([-2 -3])
