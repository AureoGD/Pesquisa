clear all
clc

% b = 2;
% a = [1 3 2];
% 
% [Ac,Bc,Cc,Dc] = tf2ss(b,a);
% sys = ss(Ac,Bc,Cc,Dc);
% sysd = c2d(sys,0.1)
% 
% A = sysd.A;
% B = sysd.B;
% C = sysd.C;

% %%
A = [0.7326 -0.0861;
     0.1722 0.9909];
B = [0.0609; 0.0064];
C = [0 1.4142];
Q = eye(2);
R = 0.01;
Ny = 2;
Nu = Ny;
Nc = 1;

%P dado pela solução de Lyapunov Eq 4
P = dlyap(A,Q)
%P = eye(2);
% H = 2*[(B'*A'*P*A*B + R + B'*Q*B)    B'*A'*P*B;
%         B'*P*A*B                    (B'*P*B+R)]
%         
% F = [((A^2)'*P*A*B + (B'*A'*P*(A^2))' + A'*Q*B + (B'*Q*A)') ((A^2)'*P*B + (B'*P*(A^2))')]

% M = [];
% phi = [];
% for i = 1:Ny
%     M = [A^(i-1)*B M];
%     phi = [phi; A^i];
% end
% omega = [B zeros(2,1) ; A*B B]
% 
% Q = eye(4);
% R = R*eye(2);
% H = 2*(M*P*M' + omega'*Q*omega + R)
% 
% F = ((A^Ny)'*P*M' + (M*P*A^(Ny))' + (omega'*Q*phi)' + (phi'*Q*omega)')

Sx = [eye(2,2) ; A ; A^2]
Su = [zeros(2,2) ;B zeros(2,1) ; A*B B]
Qlinha = blkdiag(Q,Q,P)
Rlinha = R*eye(2);

% 
% % Sx = [eye(2,2);A];
% % Su = [zeros(2,1) ; B];
% Qlinha = blkdiag(Q,P)
% Rlinha = R*eye(2)
H = (Su'*Qlinha*Su + Rlinha)
F = (Sx'*Qlinha*Su)


