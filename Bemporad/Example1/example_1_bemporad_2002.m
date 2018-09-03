clear all
close all
clc
figure(1)
hold on

%System from Example 1 Bemporad 2002
A = [0.7326 -0.0861;
     0.1722 0.9909];
B = [0.0609; 0.0064];
C = [0 1.4142];
Q = eye(2);
R = 0.01;
Ny = 2;
Nu = 2;
Nc = 2;
Umax = 2;
Umin = -2;
x_inicial = [1 1]';
%P dado pela solução de Lyapunov Eq 4
P = dlyap(A,Q);

Sx = eye(length(A));
for i = 1:Ny
    Sx = [Sx; A^i];
end

Su = zeros((Ny+1)*length(B),Nu);
for j = 1:Nu
    for i = 1:(Ny+1)
        if (i-j-1)>0
            Su((length(B)*(i-1)+1):(length(B)*(i)),j) =  [A^(i-j-1)*B];
        elseif (i-j-1) == 0
            Su((length(B)*(i-1)+1):(length(B)*(i)),j) = [B];
        end
    end
end

Qlinha = Q;
for i = 1:(Nu-1)
    Qlinha = blkdiag(Qlinha,Q);
end
Qlinha = blkdiag(Qlinha,P);
Rlinha = R*eye(Nu);

H = (Su'*Qlinha*Su + Rlinha);
F = (Sx'*Qlinha*Su);

G = [ 1  0;
     -1  0;
      0  1;
      0  -1];
W = [Umax; -Umin; Umax; -Umin];

E = zeros(Nc*Nu,Nu);

%Adição da restrição no estado 1.5 < x1,x2 < 1.5
W = [W; 1.5; 1.5; 1.5; 1.5];
E = [E; -1 0; 0 -1; 1 0; 0 1];
G = [G; zeros(4,2)];


%%Criação de CR0 e primeira divisão das areas
%% Calculo z0
S = E + G*inv(H)*F';
z0 = optimal_z_mp_QP( G, W, S, H, F, x_inicial, Nu);

%% Obtenção G_tio, W_tio e S_tio de CR0

tol = 1e-8;
[G_tio, W_tio, S_tio] = verify_active_constraints(G, W, S, x_inicial, z0, tol);

[A_CR0 b_CR0] = define_region(G, W, S, G_tio, W_tio, S_tio, H, tol);

Regions{1,1} = A_CR0;
Regions{1,2} = b_CR0;

out_X{1,1} = [1 0; 0 1; -1 0; 0 -1];
out_X{1,2} = [1.5; 1.5; 1.5; 1.5];
i=1;
Nx = 4;          %Quantidade de linhas que descrevem a restrição no espaço
%%
CR0_rest = new_rest_regions( A_CR0, b_CR0, out_X, out_X);

%%
plotregion(-A_CR0,-b_CR0)
hold on
for i = 1:3
   plotregion(-CR0_rest{i,1},-CR0_rest{i,2})
end
xlim([-1.5 1.5])
ylim([-1.5 1.5])

%%
figure(4)
hold on
new_Regions = partition_rest_regions(CR0_rest, G, W, S, H, F, tol, Nu,out_X, 1);
xlim([-1.5 1.5])
ylim([-1.5 1.5])
%%
i=1;
Rest_regions = CR0_rest;
 %for i = 1:size(Rest_regions,1)
        A = Rest_regions{i,1};
        b = Rest_regions{i,2};
        [ xc , r] = chebychev_ball( A, b, G, W, S, H, F );  %Find epsilon and x0
        plot(xc(1),xc(2),'*')
        %%
%         if r <= 1.0000e-05
%              Regions = [Regions; {Rest_regions{i,:}}];
%         else
            [ z0 ] = optimal_z_mp_QP( G, W, S, H, F, xc, Nu);  % Find z0
            [G_tio W_tio S_tio] = verify_active_constraints(G, W, S, xc, z0, tol); %Find G_tio, W_tio e S_tio
            %Define U(x)
            [ A, b ] = define_region( G, W, S, G_tio, W_tio, S_tio, H, tol) %Define polyhedron
            %%
            CR = {A b}
            Regions = [Regions; CR];
            i=1;
            %%
            nr = size(Rest_regions{i,1},1)
            Rest_CR = new_rest_regions(A,b,{Rest_regions{i,:}}, Nx+nr);
           
            %%
            Rest_regions = Rest_CR;
 %for i = 1:size(Rest_regions,1)
        A = Rest_regions{i,1};
        b = Rest_regions{i,2};
        [ xc , r] = chebychev_ball( A, b, G, W, S, H, F );  %Find epsilon and x0
        plot(xc(1),xc(2),'*')
        %%
%         if r <= 1.0000e-05
%              Regions = [Regions; {Rest_regions{i,:}}];
%         else
            [ z0 ] = optimal_z_mp_QP( G, W, S, H, F, xc, Nu);  % Find z0
            [G_tio W_tio S_tio] = verify_active_constraints(G, W, S, xc, z0, tol); %Find G_tio, W_tio e S_tio
            %Define U(x)
            [ A, b ] = define_region( G, W, S, G_tio, W_tio, S_tio, H, tol) %Define polyhedron
            %%
            CR = {A b}
            Regions = [Regions; CR];
            i=1;
            %%
            nr = size(Rest_regions{i,1},1)
            Rest_CR = new_rest_regions(A,b,{Rest_regions{i,:}}, Nx+nr);

            
            
            %%
            Nr = size(Rest_regions{i,1},1);
            partition_rest_regions(Rest_CR, G, W, S, H, F, tol, Nu, Nr);
%        end
   % end
   
   
   









