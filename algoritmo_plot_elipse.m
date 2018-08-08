clc
clear all
%Q = [1 0.5; 0.5 1];
Q = [8.78994109012211e-05 -6.476763375948896e-06; -6.476763375948896e-06 7.94625912840358e-06]

angle = 0.5*atan(2*Q(1,2)/((Q(1,1)-Q(2,2))))
% Rx = (0.5*(Q(1,1) + Q(2,2)))^0.5
% Ry = (0.5*(Q(1,1) - Q(2,2)))^0.5

% X^2/(Rx^2+Ry^2) + Y^2/(Rx^2-Ry^2)=
% 
% Rx^2+Ry^2 = Q(1,1)
% Rx^2-Ry^2 = Q(2,2)
% Rx^2 =Q(1,1)+Q(2,2)/
% 
% 2Ry^2 = Q(1,1)-Q(2,2)

m11 = Q(1,1)
m12 = Q(1,2)
m22 = Q(2,2)

Rx = (2/((m11+m22) + ((m11-m22)^2 + 4*m12^2)^0.5))^0.5
Ry = (2/((m11+m22) - ((m11-m22)^2 + 4*m12^2)^0.5))^0.5



n = 1;
for i=0:(pi/1000):(2*pi)
    x(n) = Rx*cos(i)*cos(angle)-Ry*sin(i)*sin(angle);
    y(n) = Rx*cos(i)*sin(angle)+Ry*sin(i)*cos(angle);
    n = n + 1;
end

plot(x,y,'-')
grid


%%
close all
N=1000;
Q1 = inQ_table{1,1};
Q2 = inQ_table{1,2};
Q3 = inQ_table{1,3};
Q4 = inQ_table{1,4};
Q5 = inQ_table{1,5};
plot_elipse(Q1,N);
hold on
plot_elipse(Q2,N);
%%
plot_elipse(Q3,N);
%%
plot_elipse(Q4,N);
%%
plot_elipse(Q5,N);

