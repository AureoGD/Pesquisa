clear all
close all
clc

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
P = dlyap(A,Q)

Sx = eye(length(A));
for i = 1:Ny
    Sx = [Sx; A^i];
end
Sx

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
Su

Qlinha = Q;
for i = 1:(Nu-1)
    Qlinha = blkdiag(Qlinha,Q);
end
Qlinha = blkdiag(Qlinha,P)
Rlinha = R*eye(Nu)

H = (Su'*Qlinha*Su + Rlinha)
F = (Sx'*Qlinha*Su)

G = [ 1  0;
     -1  0;
      0  1;
      0  -1]
W = [Umax; -Umin; Umax; -Umin]

E = zeros(Nc*Nu,Nu)

%Adição da restrição no estado 1.5 < x1,x2 < 1.5
W = [W; 1.5; 1.5; 1.5; 1.5];
E = [E; -1 0; 0 -1; 1 0; 0 1];
G = [G; zeros(4,2)];


%% Calculo z0
U       = sdpvar(Nu,1,'full'); %Cria o vetor das ações de controle

z = U + inv(H)*F'*x_inicial; 
S = E + G*inv(H)*F';

LMI = [];
LMI = [LMI, -G*z + W + S*x_inicial >= 0 ];

objetivo = 0.5*z'*H*z;
options = sdpsettings;
options.solver = 'sedumi';
%options.solver = 'sdpt3';

teste = solvesdp(LMI,objetivo,options);
z0 = double(z)

%% Obtenção G_tio, W_tio e S_tio de CR0

tol = 10e-6;
[G_tio, W_tio, S_tio] = verify_active_constraints(G, W, S, x_inicial, z0, tol)

[A_CR0 b_CR0] = define_region(G, W, S, G_tio, W_tio, S_tio, H)

plotregion(-A_CR0,-b_CR0)
hold on

% A_CR1 = [A_CR0(1,:); -A_CR0(2,:); A_CR0(3:length(A_CR0),:)];
% b_CR1 = [b_CR0(1,:); -b_CR0(2,:); b_CR0(3:length(b_CR0),:)];    
% plotregion(-A_CR1,-b_CR1)

Nx = 4          %Quantidade de linhas que descrevem a restrição no espaço
CR0_rest = partition_region( A_CR0, b_CR0, Nx)

for i = 1:3
    plotregion(-CR0_rest{i,1},-CR0_rest{i,2})
end

xlim([-1.5 1.5])
ylim([-1.5 1.5])


%% Chebyshev ball
xc = sdpvar(2,1);
r = sdpvar(1);
optimize(CR0_rest{1,1}*xc+r*sqrt(sum(CR0_rest{1,1}.^2,2)) <= CR0_rest{1,2},-r)
xc = double(xc);
plot(xc(1),xc(2), '*')
viscircles(xc',double(r))

xc = sdpvar(2,1);
r = sdpvar(1);
optimize(CR0_rest{2,1}*xc+r*sqrt(sum(CR0_rest{2,1}.^2,2)) <= CR0_rest{2,2},-r)
xc = double(xc);
plot(xc(1),xc(2), '*')
viscircles(xc',double(r))

xc = sdpvar(2,1);
r = sdpvar(1);
optimize(CR0_rest{3,1}*xc+r*sqrt(sum(CR0_rest{3,1}.^2,2)) <= CR0_rest{3,2},-r)
xc = double(xc)
plot(xc(1),xc(2), '*')
viscircles(xc',double(r))








