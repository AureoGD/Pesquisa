clear all
clc

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

%% Obtenção do estado factivel (x_inicial) por Chebychev Ball
% U       = sdpvar(Nu,1,'full') %Cria o vetor das ações de controle
% x0      = sdpvar(length(A),'full') %Cria o vetor dos estado
% epson   = sdpvar(1,1,'full') 
% 
% S = E + G*inv(H)*F'
% 
%% Calculo z0
U       = sdpvar(Nu,1,'full'); %Cria o vetor das ações de controle

z = U + inv(H)*F'*x_inicial; 
S = E + G*inv(H)*F';

LMI = [-G*z + W - S*x_inicial >= 0 ];

objetivo= 0.5*z'*H*z;
options=sdpsettings;
options.solver='sedumi';

teste = solvesdp(LMI,objetivo,options);
z0 = double(z)

%% Obtenção da Região Critica
G_tio = sym('G', [Ny*Nu Nu])
S_tio = sym('S', [Ny*length(A) length(A)])
W_tio = sym('W', [Ny*Nu 1])

eqn = [ G_tio*z0 - S_tio*x_inicial - W_tio == 0];
   
[HUE HUE1 HUE2 ] = solve(eqn,[G_tio S_tio W_tio])
%%


% G_tio   = sdpvar(Ny*Nu,Nu,'full');
% S_tio   = sdpvar(Ny*length(A),length(A));
% W_tio   = sdpvar(Ny*Nu,1);
% 
% EQUALITY = [ G_tio*z0 - S_tio*x0 - W_tio == 0];




%% Participação das outras regiões





%%

