function [ x_elipse, y_elipse, a, b, theta ] = plot_elipse( Q , N )
%plot_elipse Plota uma Elipse a partir da forma quadratica x'Qx = 1
%   Q - Matriz simetrica da equação
%   N - Número de pontos

    x_elipse = zeros(N);
    y_elipse = zeros(N);

    theta = 0.5*atan(2*Q(1,2)/((Q(1,1)-Q(2,2))));
    m11 = Q(1,1);
    m12 = Q(1,2);
    m22 = Q(2,2);
    a = (2/((m11+m22) + ((m11-m22)^2 + 4*m12^2)^0.5))^0.5;
    b = (2/((m11+m22) - ((m11-m22)^2 + 4*m12^2)^0.5))^0.5;

    n = 1;
    for i=0:(pi/(N/2)):(2*pi-pi/(N/2))
        x_elipse(n) = a*cos(i)*cos(theta) - b*sin(i)*sin(theta);
        y_elipse(n) = a*cos(i)*sin(theta) + b*sin(i)*cos(theta);
        n = n + 1;
    end

    plot(x_elipse,y_elipse,'-')

end

