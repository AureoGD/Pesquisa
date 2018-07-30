clc
clearvars x y x_e i ref t u e U du e e_rms
close all
Tsimu = 5;
Nint = Tsimu/Ts;

Qq = [0.1 0; 0 0.1];

Rq = [100];

Kctr = dlqr(Aad, Bad, Qq, Rq);

du(1:Nint) = 0;
U(1:Nint) = 0;
ref(1:Nint) = 350;
x(1:Nint) = 0;
dx(1:Nint) = 0;
y(1:Nint) = 0;
e(1:Nint) = 0;
w(1:Nint) = 0;
id(1:Nint)= 0;
J(1:Nint)= 0;


% for i = 1:150
% %     ref(i)= 350*i/150;
%     ref(i) = 100;
% end
% 
% for i = 1:200
% %  ref(i+150)= 350;
%    ref(i+150) = 350;
% end
%     
% for i = 1:150
% %    ref(i+350)= -350*i/150 +350;
%      ref(i+350)= -120;
% end
i = 0;

for k = 2:Nint
    %planta    
    x(k)   = sysd.A*x(k-1) + sysd.B*U(k-1); 
    y(k)   = sysd.C*x(k);
    
    %obtenção da derivada do estado, do erro, e montagem do vetor x aumentado
    dx(k) = x(k) - x(k-1);
    e(k)  = y(k) - ref(k);
    x_a    = [dx(k); e(k)];% este sinal deveria ser positivo
    
    U(k) = U(k-1)-Kctr*x_a;
end

t = 0:(Nint-1);
t = t*Ts;
plot(t,y(1,:))
hold on
plot(t,ref(1,:))
title('Output')

figure()
plot(t,U(1,:))
title('U')

figure()
plot(t,du(1,:))
title('Delta U')

figure()
e_rms(1:Nint)= 0;
for i=1:Nint
    e_rms(i) = sqrt(e(i)^2); 
end
plot(t,e_rms(1,:))

