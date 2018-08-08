clc
close all
clear all

%Sistema Discreto:
A = [1 0.1;
    0 1];
B = [1 ;
    0.1*0.787];
C = [1 0];

Ts = 0.1;

A1 = [1 0.1
      0  0.99];
A2 = [1 0.1
      0   0];

%Parametros do Controlador
U_MAX = 2;
Y_MAX = 0.2;
R = 0.00002;
Qq = [1 0; 
      0 0];
% Número de Elipses
S = 13;
%Estados Factiveis do Sistema
X1 = [1, 0.1574;
      0.9, 0.9*0.1574;
      0.75, 0.75*0.1574;
      0.65, 0.65*0.1574;
      0.52, 0.52*0.1574;
      0.4,  0.4*0.1574; 
      0.28, 0.28*0.1574;
      0.18, 0.18*0.1574;
      0.1,  0.1*0.1574;
      0.05, 0.05*0.1574;
      0.02, 0.02*0.1574;
      0.01, 0.01*0.1574;
      0.001, 0.001*0.1574];
X1(:,2)=2*0.1574;
 X1 = X1';
 %% Criação das LMI e execucao do algortimo offline 
clc
clearvars gamma_table Q_table inQ_table F_table Ur_table Y_table Z_table
 
  

for i = 1:S
    clearvars gamma Q F Y Z Xrr 
    
    %Variaveis de busca:
    gamma = sdpvar(1,1,'full');
    Q     = sdpvar(2,2,'symmetric'); %Define automaticamente que é simetrica
    % F obtida por F = Y*inv(Q);
    Xrr   = sdpvar(1,1,'full');
    Y     = sdpvar(1,2,'full'); %Define assim por padrão
    %Z   = sdpvar(1,1,'full');

    LMI1 = [1 , X1(:,i)';
            X1(:,i), Q];
    
    LMI2_1 = [    Q,      (Q*(A1')+(Y')*(B')),      Q*(Qq^0.5),       (Y')*(R^0.5);
                A1*Q+B*Y,         Q,           zeros(2,2),        zeros(2,1);
               (Qq^0.5)*Q,   zeros(2,1),      gamma*eye(2),      zeros(2,2);
               (R^0.5)*Y,    zeros(1,2),      zeros(1,2),      gamma*eye(1)];
     
    LMI2_2 = [    Q,      (Q*(A2')+(Y')*(B')),      Q*(Qq^0.5),       (Y')*(R^0.5);
                A2*Q+B*Y,         Q,           zeros(2,2),        zeros(2,1);
               (Qq^0.5)*Q,   zeros(2,1),      gamma*eye(2),      zeros(2,2);
               (R^0.5)*Y,    zeros(1,2),      zeros(1,2),      gamma*eye(1)];
    
    Ur = U_MAX^2;
    LMI3 = [  Xrr,     Y;
              Y',     Q];
       
    LMI5 = Ur - Xrr;
      
%     Z = Y_MAX^2;    
%     LMI4_1 = [       Z,           C*(A1*Q + B*Y);
%              ((A1*Q + B*Y)')*(C'),           Q    ];
%     
%     LMI4_2 = [       Z,           C*(A2*Q + B*Y);
%              ((A2*Q + B*Y)')*(C'),           Q    ];   
%         

    LMIs = [];
    LMIs = [LMIs,    Q    >=0];
    LMIs = [LMIs,  LMI1   >=0];
    LMIs = [LMIs,  LMI2_1 >=0];
    LMIs = [LMIs,  LMI2_2 >=0];
    LMIs = [LMIs,  LMI3   >=0];
    LMIs = [LMIs,  LMI5   >=0];
%     LMIs = [LMIs,  LMI4_1 >=0];
%     LMIs = [LMIs,  LMI4_2 >=0];
    if (i>1)
        LMIs  =  [LMIs , (Q_table{i-1} - Q) >= 0];
        %LMIs  =  [LMIs , X1(i,:)*Q_table{i-1}*X1(i,:)' - X1(i,:)*Q*X1(i,:)' > 0];
    end
    
    objetivo=[gamma];
    options=sdpsettings;
    option.solver='sedumi';
    solvesdp(LMIs,objetivo,options);
   
    Q = double(Q);
    Z = 1;%double(Z);
    Y = double(Y);
    invQ = inv(Q);
    F = Y*invQ;
    Xrr = double(Xrr);
    
    gamma_table{i} = gamma;
    Q_table{i} = Q;
    inQ_table{i} = invQ;
    F_table{i} = F;
    Ur_table{i} = Ur;
    Y_table{i} = Y;
    Z_table{i} = Z;  
    Xrr_table{i} = Xrr;
end