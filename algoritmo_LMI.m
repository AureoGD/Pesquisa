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

%Sistema aumenado:
Aad = [sysd.A, 0; Ts*sysd.C 1];
Bad = [sysd.B; 0];

%limitação maxima do sistema:
y_max = 350;
x_max = y_max/C;
u_max  = 100;


%% LMI's para a solução do off line-MPC: Função custo J = x'Qx+u'Ru com estado aumentado
clc
clearvars X Y_table Q_table inQ_table F_table
Aad1 = [sysd.A*0.05 0;
        Ts*sysd.C   1]; %criação do politopo convexo de Aad2
    
Aad2 = [sysd.A*1.05 0;
        Ts*sysd.C   1]; %criação do politopo convexo de Aad2

% Cad = [1 1];
% %Ponderação de estados e ação de controle:
% Y_max = [0.0944 ; C*x_max].^2;

%U_max = [u_max]^2;

U_max = [10]^2;

% É preciso encontra S estados do sistema, e criar um vetor com eles
S = 5 %nº de elipses

X = [x_max   , C*x_max;
     x_max/2 , C*x_max/2;
     x_max/4 , C*x_max/4;
     x_max/16, C*x_max/26;
     x_max/32, C*x_max/32];
 
%Qq    = [1/(sysd.B*u_max)^2 0; 
 %        0 1/(350)^2];
     
%Rq    = [1/u_max^2];

Qq = [1/x_max 0; 0 1/C*x_max];

Rq = [1];

for i=1:S
    
    clearvars Y Q U_r Y_r gamma F
    
    %Variaveis de busca:
    Y     = sdpvar(1,2,'full'); %Define assim por padrão
    Q     = sdpvar(2,2,'symmetric'); %Define automaticamente que é simetrica
    U_r   = sdpvar(1,1,'full');
    Y_r   = sdpvar(1,1,'full');
    gamma = sdpvar(1,1,'full');
  
    LMI1   = [  1      ,  X(i,:)
               X(i,:)' ,   Q ];
    
    LMI2_1 = [     Q        , (Q*Aad1'+Y'*Bad') , Q*Qq^1/2     , Y'*Rq^1/2;
             (Aad1*Q+Bad*Y) ,         Q         , zeros(2,2)   , zeros(2,1);      
              Qq^1/2*Q      ,    zeros(2,1)     , gamma*eye(2) , zeros(2,2);
              Rq^1/2*Y      ,    zeros(1,2)     , zeros(1,2)   , gamma*eye(1)];
   
    LMI2_2 = [     Q        , (Q*Aad2'+Y'*Bad') , Q*Qq^1/2     , Y'*Rq^1/2;
             (Aad2*Q+Bad*Y) ,         Q         , zeros(2,2)   , zeros(2,1);      
              Qq^1/2*Q      ,    zeros(2,1)     , gamma*eye(2) , zeros(2,2);
              Rq^1/2*Y      ,    zeros(1,2)     , zeros(1,2)   , gamma*eye(1)];
           
     Res_U = [U_r, Y;
              Y' , Q];
          
%      Res_Y1 = [     Q              ,   (Aad1*Q+ Bad*Y)'*Cad';
%                Cad*(Aad1*Q+ Bad*Y) ,           Y_r        ];
%           
%      Res_Y2 = [     Q              ,   (Aad2*Q+ Bad*Y)'*Cad';
%                Cad*(Aad2*Q+ Bad*Y) ,           Y_r        ];
    
    %verificar o help do sdpvar
        LMIs  =  [];
        LMIs  =  [LMIs , Q           > 0];
        LMIs  =  [LMIs , LMI1        >= 0];
        LMIs  =  [LMIs , LMI2_1      >= 0];
        LMIs  =  [LMIs , LMI2_2      >= 0];
        LMIs  =  [LMIs , Res_U       >= 0];
        LMIs  =  [LMIs , U_max-U_r   >= 0];
%         LMIs  =  [LMIs , Res_Y1    >= 0];
%         LMIs  =  [LMIs , Res_Y2    >= 0];
%         LMIs  =  [LMIs , Y_max-Y_r >= 0];
    if (i~=1)
        LMIs  =  [LMIs , Q_table{i-1} - Q > 0];
    end
    objetivo=[gamma];
   
    options=sdpsettings;
    option.solver='sedumi';

    solvesdp(LMIs,objetivo,options);

    Q = double(Q)
    Y = double(Y)

    F = Y*inv(Q);
    
     Y_table{i}   = Y;
     F_table{i}   = F;
     Q_table{i}   = Q;
     inQ_table{i} = inv(Q);
end

