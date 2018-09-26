%SetPoint Control
%System from Example 1 Bemporad 2002
Ap = [0.7326 -0.0861;
     0.1722 0.9909];
Bp = [0.0609; 0.0064];
Cp = [0 1.4142];

A = [Ap zeros(1,2)';
     Cp*Ap eye(1)]
B = [Bp;
    Cp*Bp]
C = [zeros(1,2) eye(1)];

%A = [A zeros(2,1); C 0];
%B = [B ; 0];
Nint = 70;

x1 = zeros(1,Nint);
x2 = x1;
y = zeros(Nint,1);
u = zeros(Nint,1);
delta_u = zeros(Nint,1);
epson = zeros(Nint,1);
r = 0.5;
x0 = [1 1];
X(:,1) = [ x0 Cp*x0'];
X(:,2) = [ x0 Cp*x0'];
x1(1) = x0(1);
x1(2) = x0(1);
x2(1) = x0(2);
x2(2) = x0(2);
epson(1) = 0;

%index = 1;
[X(:,2) ; u(1) ; r]
for i = 2:Nint
    for j = 1:size(Regions,1)
        A_CRi = Regions{j,1};
        b_CRi = Regions{j,2};
        flag = 0;
        for k = 1:size(A_CRi,1)
            if(A_CRi(k,:)*[X(:,i) ; u(i-1) ; r] > b_CRi(k))
                flag = 1;
            end
        end
        if flag == 0
            index = j
        end
    end
    indices(i) = index;
    delta_u(i) = Regions{index,3}*[X(:,i) ; u(i-1) ; r] + Regions{index,4};
%     A*x(:,i)
%     B*u
    
%    %DU = quadprog(H,([X(:,i) ; u(i-1) ; r]'*F)');
%     DU = sdpvar(Nu,1,'full');
%     z = DU + inv(H)*F'*[X(:,i) ; u(i-1) ; r]; 
% 
%     LMI = [];
%     %LMI = [LMI, -G*z + W + S*x0 >= 0 ];
%     LMI = [LMI, G*z <= W + S*[X(:,i) ; u(i-1) ; r]];
%     objetivo = 0.5*z'*H*z;
%     options = sdpsettings;
%     options.solver = 'sedumi';
%     options.verbose = 0;
%     diagnostics = optimize(LMI,objetivo,options);
%     DU_hist(:,i) = double(DU);
    
    %delta_u(i) = double(DU(1));


    X(:,i+1) = A*X(:,i)+B*delta_u(i);
    y(i) = C*X(:,i);
    
    u(i) = delta_u(i) + u(i-1);
    y(i+1) = X(3,i);
    x1(i+1) = X(1,i);
    x2(i+1) = X(2,i);
    %epson(i) = r - y(i);

end

%%
figure(10)
stairs(x1)
hold on
stairs(x2,'r')
legend('x1','x2')

figure(11)
stairs(u)
legend('u')

figure(12)
stairs(y)
legend('y')

figure(13)
stairs(indices)
legend('index')

figure(14)
stairs(delta_u)
legend('Delta U')

