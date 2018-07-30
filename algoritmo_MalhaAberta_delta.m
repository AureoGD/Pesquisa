clc
clearvars x y x_e i ref t u e dx
close all
Tsimu = 1;
Nint = Tsimu/Ts;

du(1:Nint) = 0;
U(1:Nint) = 100;
ref(1:Nint) = 350;
x(1:Nint) = 0;
dx(1:Nint) = 0;
y(1:Nint) = 0;
e(1:Nint) = 0;
w(1:Nint) = 0;
id(1:Nint)= 0;

% for i = 1:150
%     ref(i)= 350*i/150;
% end
% 
% for i = 1:200
%    ref(i+150)= 350; 
% end
%     
% for i = 1:150
%     ref(i+350)= -350*i/150 +350;
% end
i = 0;

for k = 2:Nint
    %planta    
    x(k)   = sysd.A*x(:,k-1)+ sysd.B*U(k-1); 
    dx(k) = x(k)- x(k-1); 
end

t = 0:(Nint-1);
t = t*Ts;
plot(t,dx(1,:))
title('xd state')
figure()
plot(t,x(1,:))
title('Output state')