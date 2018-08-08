clc
close all
clear all

%Definição do sistema continuo:
A = [-11.58];
B = [1];
C = [42.74];
D = [0];

sys = ss(A, B, C, D);

%Definição do sistema discreto:

Ts = 0.01;
sysd = c2d(sys,Ts);

%Sistema aumentado:
% Aad = [sysd.A,         0; 
%       Ts*sysd.C*sysd.A 1];
% 
% Bad = [sysd.B; 
%        sysd.C*sysd.B];

Ad = eye(size(A))+A*Ts
Bd = B*Ts
Cd = C
   
%Funcao que calcula as matrizes do modelo aumentado e todas as matrizes do
%MPC - mpcgain(Ap,Bp,Cp,Np,Nc);

[m1,n1]=size(Cd);
[n1,n_in]=size(Bd);
A_e=eye(n1+m1,n1+m1);
A_e(1:n1,1:n1)=Ad;
A_e(n1+1:n1+m1,1:n1)=Cd*Ad;
B_e=zeros(n1+m1,n_in);
B_e(1:n1,:)=Bd;
B_e(n1+1:n1+m1,:)=Cd*Bd;
C_e=zeros(m1,n1+m1);
C_e(:,n1+1:n1+m1)=eye(m1,m1);


%%
%limitação maxima do sistema:
y_max = 350;
y_min = -100;
x_max = y_max/C;
Umax  = 200;
Umin = -50;
delta_Umax = 10;
delta_Umin = -10;

%% Inicizalição de vetores
Tsimu = 10;
Nint = Tsimu/Ts;

du(1:Nint) = 0;
U(1:Nint) = 0;
x(1:Nint) = 0;
dx(1:Nint) = 0;
y(1:Nint) = 0;
e(1:Nint) = 0;
w(1:Nint) = 0;
id(1:Nint)= 0;
J(1:Nint)= 0;
x_a(2,1:Nint) = 0;
Xf = [0 0]';

%Criação referencia
ref(1:Nint)= 10;
% for i = 1:300
% %     ref(i)= 350*i/150;
%     ref(i) = 100;
% end
% 
% for i = 1:400
% %  ref(i+150)= 350;
%    ref(i+300) = 350;
% end
%     
% for i = 1:300
% %    ref(i+350)= -350*i/150 +350;
%      ref(i+700)= -120;
% end

%% Paramentros MPC e ponderações
Np = 20;
Nc = 15;
Qq = eye(2);
Rq = [(1/100)^2];
R_bar = Rq*eye(Nc);

%% Calculo Matrizes MPC
[m1,n1]=size(C_e);
n=n1+m1;
h(1,:)=C_e;
F(1,:)=C_e*A_e;
for kk=2:Np
h(kk,:)=h(kk-1,:)*A_e;
F(kk,:)= F(kk-1,:)*A_e;
end
v=h*B_e;
Phi=zeros(Np,Nc); %declare the dimension of Phi
Phi(:,1)=v; % first column of Phi
for i=2:Nc
Phi(:,i)=[zeros(i-1,1);v(1:Np-i+1,1)]; %Toeplitz matrix
end
BarRs=ones(Np,1);
Phi_Phi= Phi'*Phi;
Phi_F= Phi'*F;
Phi_R=Phi'*BarRs;

% Restrições
Ac = zeros(4*Nc-2,Nc);
b0 = zeros(4*Nc-2,1);
% for k = 1:Nc
%     Ac(4*k-3,k) = 1;
%     b0(4*k-3,1) = -delta_Umin;
%     Ac(4*k-2,k) = -1;
%     b0(4*k-2,1)= delta_Umax;
%     if k < Nc
%         Ac(4*k-1,k) = -1;
%         Ac(4*k-1,k+1) = 1;
%         b0(4*k-1,1) = -delta_Umin;
%         Ac(4*k,k) = 1;
%         Ac(4*k,k+1) = -1;
%         b0(4*k,1) = - delta_Umax;
%     end
% end
% Ac = Ac
% b0 = b0
triang_Nc = tril(ones(Nc));
Ac = [   eye(Nc);
        -eye(Nc);
       triang_Nc;
      -triang_Nc;
             Phi;
            -Phi;
      ]
b0 = [ones(Nc,1)*delta_Umax;
      -ones(Nc,1)*delta_Umin;
      ones(Nc,1)*delta_Umax;
      -ones(Nc,1)*delta_Umin;
      ones(Np,1)*y_max - F*Xf;
      F*Xf - ones(Np,1)*y_min;
];

%% Simulação

du_Nc(1:Nc) = 0;
for k=2:Nint
%     x(k)   = sysd.A*x(k-1) + sysd.B*U(k-1); 
%     y(k)   = sysd.C*x(k);
    x(k)   = Ad*x(k-1) + Bd*U(k-1); 
    y(k)   = Cd*x(k);
    dx(k) = x(k) - x(k-1);
    e(k)  = y(k) - ref(k);
    
    %Xf = [dx(k) y(k)]';
    Xf = [x(k) U(k-1)]';
    
    %Atualiza restrição
%     b0(1,1) = U(k-1) + delta_Umin;
%     b0(2,1) = -U(k-1) - delta_Umax;
b0 = [ones(Nc,1)*delta_Umax;
      -ones(Nc,1)*delta_Umin;
      ones(Nc,1)*(delta_Umax-U(k-1));
      -ones(Nc,1)*(delta_Umin-U(k-1));
      ones(Np,1)*y_max - F*Xf;
      F*Xf - ones(Np,1)*y_min;
];


    x_a(:,k) = [dx(k); -e(k)];
    
    Hqp = double(Phi_Phi+R_bar);
    Fqp = double(Phi_F*Xf-Phi_R*ref(k));
   
%     %Solucao usando mpcqpsolver
    [L, p] = chol(Hqp,'lower');
    Linv = L\eye(size(Hqp,1));

    opt_mpc = mpcqpsolverOptions;
    opt_mpc.IntegrityChecks=false;
    opt_mpc.MaxIter=10;
    opt_mpc.UseSuboptimalSolution = true;

    iA = false(size(b0));
    %[du(k), status, iA] = mpcqpsolver(Linv,Fqp,Ac,b0,[], zeros(0,1),iA,opt_mpc);
    %[du_Nc, status, iA] = mpcqpsolver(Linv,Fqp,Ac,b0,[],zeros(0,1),iA,opt_mpc);
    du_Nc = quadprog(Hqp,Fqp,Ac,b0)
    du(k) = du_Nc(1,1);
    %du(k) = 0;
    %U(k) = U(k-1) + du(k);
    U(k) = du(k);
    k
end

figure(1)
t = 0:(Nint-1);
t = t*Ts;
plot(t,y(1,:))
hold on
plot(t,ref(1,:))
title('Output')

figure(2)
plot(t,U(1,:))
title('U')

figure(3)
plot(t,du(1,:))
title('Delta U')
