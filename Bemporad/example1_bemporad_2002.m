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

LMI = [];
LMI = [LMI, -G*z + W + S*x_inicial >= 0 ];

objetivo = 0.5*z'*H*z;
options = sdpsettings;
options.solver = 'sedumi';
%options.solver = 'sdpt3';

teste = solvesdp(LMI,objetivo,options);
z0 = double(z)
lambda0 = (-inv(H)*G')\z0

%%
% H_example_double = [0.8365 0.3603; 0.3603 0.2059];
% G_example_double = [1 0 ; -1  0; 0 1; 0 -1];
% F_example_double = [0.4624 1.2852; 0.1682 0.5285];
% S_example_double = zeros(4,2);%%G_example_double*inv(H_example_double)*F_example_double';
% 
% H = H_example_double;
% G = G_example_double;
% F = F_example_double;
% S = S_example_double;

%% Obtenção de G_tio, W_tio e S_tio

tol = 10e-6;
index = [];
for i = 1:(Nc*Nu)
    if ((G(i,:)*z0 - W(i,:) - S(i,:)*x_inicial < tol) && (double(G(i,:)*z0 - W(i,:) - S(i,:)*x_inicial)> -tol))
        index = [index , i];
        G(i,:)*z0 - W(i,:) - S(i,:)*x_inicial
    end
end
index
G_tio = [];
S_tio = [];
W_tio = [];

if (length(index)>0)
    for i = 1:length(index)
        G_tio = [G_tio; G(index(i),:)];
        S_tio = [S_tio; S(index(i),:)];
        W_tio = [W_tio; W(index(i),:)];
    end
end

T = (inv(H)*G_tio'*inv(G_tio*inv(H)*G_tio'));

A0_1 = G*T*S_tio-S
b_1 = W - G*T*W_tio
A0_2 = inv(G_tio*inv(H)*G_tio')*S_tio
b_2 = -inv(G_tio*inv(H)*G_tio')*W_tio

A_plot = [A0_1(3:4,:); A0_2];
b_plot = [b_1(3:4,:); b_2];

figure(1)
plotregion(-A_plot,-b_plot)

%%
figure(4)
%A_plot_2 = [3.4155 -4.6452; -0.1044 -0.1215; -0.1559 -0.0922; 1 0; 0 1]%; 1 1];
%b_plot_2 = [2.6341; -0.0353;-0.0267; 1.5; 1.5 ]%; -3];
A_plot_2 = [A_plot;  1 0; 0 1];
B_plot_2 = [b_plot; 1.5;1.5];
hold on
grid

plotregion(-A_plot_2,-b_plot_2)


%figure(3)
A_plot_R1 = [-5.922 -6.8883; 5.9220 6.8883; -1.5379 6.8291; 1.5379 -6.8291];
b_plot_R1 = 2*ones(4,1);
plotregion(-A_plot_R1,-b_plot_R1)
xlim([-1.5 1.5])
ylim([-1.5 1.5])

figure(2)
A_plot_R24 = [-3.4155  4.6452;  0.1044 0.1215; 0.1559 0.0922];
b_plot_R24 = [2.6341; -0.0353; -0.0267];
plotregion(-A_plot_R24,-b_plot_R24)
hold on
A_plot_R6 = [-6.4159 -4.6953; -0.275 0.1220; 6.4159 4.6953];
b_plot_R6 = [1.3577; -0.0357; 2.6423];
plotregion(-A_plot_R6,-b_plot_R6)

plotregion(A_plot_R6,-b_plot_R6)
xlim([-1.5 1.5])
ylim([-1.5 1.5])

% figure(5)
% V = con2vert(A_plot_R24,b_plot_R24)
% plot(V(1,:),V(2,:))

% 
% lambda_tio = -inv(G_tio*inv(H)*G_tio')*(W_tio + S_tio*x_inicial)
% 
% ganho_U = inv(H)*G_tio'*inv(G_tio*inv(H)*G_tio')*W_tio
% ganho_U = ganho_U(1,:)
% const_U = (inv(H)*G_tio'*inv(G_tio*inv(H)*G_tio')*S_tio - inv(H)*F')
% const_U = const_U(1,:)
% 
% T = (inv(H)*G_tio');
% A_CR0 = (G*T*inv(G_tio*T)*S_tio-S)
% b_CR0 = W-G*T*inv(G_tio*T)*W_tio
% 
% A_CR0_2 = inv(G_tio*T)*S_tio
% b_CR0_2 = -inv(G_tio*T)*W_tio


%% Obtenção da Região Critica
% %G_tio = sym('G', [Ny*Nu Nu])
% S_tio = sym('S', [Ny*length(A) length(A)])
% W_tio = sym('W', [Ny*Nu 1])
% solG = zeros(Ny*Nu,Nu);
% solS = zeros(Ny*length(A),length(A));
% solW = zeros(Ny*Nu,1); 
% eq1 = [G*z0 - S_tio*x_inicial - W_tio == 0 ];
% %   eq2 = [inv(H)*G_tio'*inv(G_tio*inv(H)*G_tio')*(W_tio + S_tio*x_inicial)]    ; 
% hue = solve(eq1,[S_tio W_tio],'ReturnConditions', true)
% 

% %%
% % G_tio = sdpvar(Ny*Nu,,Nu,'full')
% % S_tio = sdpvar(Ny*length(A), length(A), 'full')
% % W_tio = sdpvar(Ny*Nu,1,'full');
% % 
% % LMIs = [];
% % LMIs = [LMIs, G_tio*z0 - S_tio*x_inicial - W_tio == 0]
% % LMIs = [LMIs, W + S*x_inicial - G*inv(H)*G_tio'*inv(G_tio*inv(H)*G_tio')*(W_tio + S_tio*x_inicial) >= 0]
% % LMIs = [LMIs, -inv(G_tio*inv(H)*G_tio')*(W_tio + S_tio*x_inicial) >= 0]
% 
% %%
% 
% 
% % G_tio   = sdpvar(Ny*Nu,Nu,'full');
% % S_tio   = sdpvar(Ny*length(A),length(A));
% % W_tio   = sdpvar(Ny*Nu,1);
% % 
% % EQUALITY = [ G_tio*z0 - S_tio*x0 - W_tio == 0];
% 
% 
% 
% 
% %% Participação das outras regiões