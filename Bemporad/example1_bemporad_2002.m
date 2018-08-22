clear all
close all
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

%Adição da restrição no estado 1.5 < x1,x2 < 1.5
W = [W; 1.5; 1.5; 1.5; 1.5];
E = [E; -1 0; 0 -1; 1 0; 0 1];
G = [G; zeros(4,2)];

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

%% Obtenção de G_tio, W_tio e S_tio, e CR0

tol = 10e-6;
index = [];
for i = 1:length(W)%(Nc*Nu)
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

%Inequações que descrevem a região CR0 pelas inequações que garantem a
%resposta ótima
A_CR0_1 = G*T*S_tio-S
b_CR0_1 = W - G*T*W_tio
A_CR0_2 = inv(G_tio*inv(H)*G_tio')*S_tio
b_CR0_2 = -inv(G_tio*inv(H)*G_tio')*W_tio

%Restrição do tamanho dos estado -1.5 <= x1,x2 <= 1.5
% A_espaco = [1 0; 0 1; -1 0; 0 -1];
% b_espaco = [1.5; 1.5; 1.5; 1.5];

A_CR0 = [];
b_CR0 = [];
for i = 1:size(A_CR0_1,1)
    %if (A_CR0_1(i,:) ~= zeros(1,length(A)))
    if ((A_CR0_1(i,1) ~= 0) || (A_CR0_1(i,2) ~= 0))
        A_CR0 = [A_CR0; A_CR0_1(i,:)];
        b_CR0 = [b_CR0; b_CR0_1(i,:)];
    end
    i
end
for i = 1:size(A_CR0_2,1)
    %if (A_CR0_2(i,:) ~= zeros(1,length(A)))
    if ((A_CR0_2(i,1) ~= 0) || (A_CR0_2(i,2) ~= 0))
        A_CR0 = [A_CR0; A_CR0_2(i,:)];
        b_CR0 = [b_CR0; b_CR0_2(i,:)];
    end
end
%Adição da restrição do espaço
% A_CR0 = [A_CR0 ; A_espaco]
% b_CR0 = [b_CR0 ; b_espaco]

A_CR1 = [-A_CR0(1,:); A_CR0(3:6,:)]%;A_CR0((size(A_CR0,1)-3):size(A_CR0,1),:)];
b_CR1 = [-b_CR0(1,:); b_CR0(3:6,:)]%;b_CR0((size(b_CR0,1)-3):size(b_CR0,1),:)];

A_CR2 = [A_CR0(1,:); -A_CR0(2,:); A_CR0(3:6,:); ]%A_CR0((size(A_CR0,1)-3):size(A_CR0,1),:)];
b_CR2 = [b_CR0(1,:); -b_CR0(2,:); b_CR0(3:6,:); ]%b_CR0((size(b_CR0,1)-3):size(b_CR0,1),:)];

A_CR3 = [A_CR0(1:6,:); -A_CR0(7,:)];% A_CR0(4:7,:); ]%A_CR0((size(A_CR0,1)-3):size(A_CR0,1),:)];
b_CR3 = [b_CR0(1:6,:); -b_CR0(7,:)];% b_CR0(4:7,:); ]%b_CR0((size(b_CR0,1)-3):size(b_CR0,1),:)];

figure(1)
plotregion(-A_CR0,-b_CR0)
hold on
xlim([-1.5 1.5])
ylim([-1.5 1.5])
%%
plotregion(-A_CR1,-b_CR1)
%%
plotregion(-A_CR2,-b_CR2)
%%
plotregion(-A_CR3,-b_CR3)
grid
xlim([-1.5 1.5])
ylim([-1.5 1.5])


% --------------- Até aqui ok ----------------------------------------

%% Partição das novas regiões

z_CR1  = sdpvar(Nu,1,'full');
x_CR1  = sdpvar(length(A),1,'full');
epsilon_CR1 = sdpvar(1,1,'full');

LMI = [];
LMI = [LMI, (-G*z_CR1 + W + S*x_CR1) >= 0 ];
LMI = [LMI, (-A_CR1(1,:)*x_CR1 - epsilon_CR1*norm(A_CR1(1,:)) + b_CR1(1,:)) >= 0];

objetivo = -epsilon_CR1;
options = sdpsettings;
options.solver = 'sedumi';

teste_CR1 = solvesdp(LMI,objetivo,options);

z0_CR1 = double(z_CR1);
x0_CR1 = double(x_CR1);

index_R1 = [];
for i = 1:(Nc*Nu)
    if ((G(i,:)*z0_CR1 - W(i,:) - S(i,:)*x0_CR1 < tol) && (double(G(i,:)*z0_CR1 - W(i,:) - S(i,:)*x0_CR1)> -tol))
        index_R1 = [index_R1 , i];
    end
end
G_tio_CR1 = [];
S_tio_CR1 = [];
W_tio_CR1 = [];

if (length(index_R1)>0)
    for i = 1:length(index_R1)
        G_tio_CR1 = [G_tio; G(index_R1(i),:)];
        S_tio_CR1 = [S_tio; S(index_R1(i),:)];
        W_tio_CR1 = [W_tio; W(index_R1(i),:)];
    end
end

T = (inv(H)*G_tio_CR1'*inv(G_tio_CR1*inv(H)*G_tio_CR1'));

%Inequações que descrição a região CR0 pelas inequações que garantem a
%resposta ótima
A_CR1_1 = G*T*S_tio_CR1-S
b_CR1_1 = W - G*T*W_tio_CR1
A_CR1_2 = inv(G_tio_CR1*inv(H)*G_tio_CR1')*S_tio_CR1
b_CR1_2 = -inv(G_tio_CR1*inv(H)*G_tio_CR1')*W_tio_CR1

A_CR1_NOVO = [];
b_CR1_NOVO = [];
for i = 1:size(A_CR0_1,1)
    if (A_CR1_1(i,:) ~= zeros(1,length(A)))
        A_CR1_NOVO = [A_CR1_NOVO; A_CR1_1(i,:)];
        b_CR1_NOVO = [b_CR1_NOVO; b_CR1_1(i,:)];
    end
end
for i = 1:size(A_CR1_2,1)
    if (A_CR1_2(i,:) ~= zeros(1,length(A)))
        A_CR1_NOVO = [A_CR1_NOVO; A_CR1_2(i,:)];
        b_CR1_NOVO = [b_CR1_NOVO; b_CR1_2(i,:)];
    end
end


%%
A_CR1_NOVO = [A_CR1_NOVO(3:4,:); A_espaco];
b_CR1_NOVO = [b_CR1_NOVO(3:4,:); b_espaco];

figure(2)
plotregion(-A_CR1_NOVO,-b_CR1_NOVO)



%%
z_sao_thomeh = sdpvar(Nu,1,'full');

x0_xunxo = [0.5; 2]*10^5;
LMI = [];
LMI = [LMI, -G*z_sao_thomeh + W + S*x0_xunxo >= 0 ];

objetivo = 0.5*z_sao_thomeh'*H*z_sao_thomeh;
options = sdpsettings;
options.solver = 'sedumi';
%options.solver = 'sdpt3';

teste_sao_thomeh = solvesdp(LMI,objetivo,options);

U_sao_thomeh = double(z_sao_thomeh) - inv(H)*F'*x0_xunxo  
U_CR1 = double(z0_CR1 - inv(H)*F'*x0_xunxo)

figure(5)
viscircles(x0_xunxo', double(epsilon_CR1))
hold on
A_XUNXO_espaco = [1 0; 0 1; -1 0; 0 -1];
b_XUNXO_espaco = [3*10^5; 3*10^5; 3*10^5; 3*10^5];
A_CR1_XUNXO = [A_CR1_NOVO(1:2,:); A_XUNXO_espaco ];
b_CR1_XUNXO = [b_CR1_NOVO(1:2,:); b_XUNXO_espaco ];
plotregion(-A_CR1_XUNXO,-b_CR1_XUNXO)
xlim([-Inf Inf])
ylim([-Inf Inf])
grid

% 
% figure(1)
% A_plot = [A_CR0_1(3:4,:); A_CR0_2];
% b_plot = [b_CR0_1(3:4,:); b_CR0_2];
% %A_plot_2 = [3.4155 -4.6452; -0.1044 -0.1215; -0.1559 -0.0922; 1 0; 0 1]%; 1 1];
% %b_plot_2 = [2.6341; -0.0353;-0.0267; 1.5; 1.5 ]%; -3];
% A_plot_2 = [A_plot;  1 0; 0 1; -1 0; 0 -1];
% b_plot_2 = [b_plot; 1.5;1.5; 1.5 ; 1.5];
% grid
% plotregion(-A_plot,-b_plot)
% hold on
% 
% plotregion(-A_plot_2,-b_plot_2)
% 
% A_plot_R1 = [-5.922 -6.8883; 5.9220 6.8883; -1.5379 6.8291; 1.5379 -6.8291];
% b_plot_R1 = 2*ones(4,1);
% plotregion(-A_plot_R1,-b_plot_R1)
% 
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])







% %figure(3)
% A_plot_R1 = [-5.922 -6.8883; 5.9220 6.8883; -1.5379 6.8291; 1.5379 -6.8291];
% b_plot_R1 = 2*ones(4,1);
% plotregion(-A_plot_R1,-b_plot_R1)
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])
% 
% figure(2)
% A_plot_R24 = [-3.4155  4.6452;  0.1044 0.1215; 0.1559 0.0922];
% b_plot_R24 = [2.6341; -0.0353; -0.0267];
% plotregion(-A_plot_R24,-b_plot_R24)
% hold on
% A_plot_R6 = [-6.4159 -4.6953; -0.275 0.1220; 6.4159 4.6953];
% b_plot_R6 = [1.3577; -0.0357; 2.6423];
% plotregion(-A_plot_R6,-b_plot_R6)
% 
% plotregion(A_plot_R6,-b_plot_R6)
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])

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


%% Participação das outras regiões


