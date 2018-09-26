clear all
close all
clc

Ap = [0.7326 -0.0861;
     0.1722 0.9909];
Bp = [0.0609; 0.0064];
Cp = [0 1.4142];

A = [Ap zeros(1,2)';
     Cp*Ap eye(1)]
B = [Bp;
    Cp*Bp]
C = [zeros(1,2) eye(1)];

Ny = 3;
Nu = 3;
Nc = 3;

Q = eye(1);
R = 0.01;

Rlinha = [];
for i = 1:Nu
    Rlinha = blkdiag(Rlinha,R);
end

Qlinha = [];
for i = 1:Ny
    Qlinha = blkdiag(Qlinha,Q);
end

model = LTISystem('A',A,'B',B,'C',C);

%%
model.x.min = [-10; -10; -15];
model.x.max = [10; 10; 15];
model.u.min = [-50];
model.u.max = [50];
%%
model.x.with('reference');
model.x.reference = 'free';
model.x.penalty = QuadFunction(Qlinha);

%%
model.u.with('deltaPenalty');
model.u.deltaPenalty = QuadFunction(R);

ctrl = MPCController(model, Ny);
%%
ectrl = ctrl.toExplicit();

%%
loop = ClosedLoop( ectrl, model);
x0 = [1 1 1.44]';
u0 = 0;
Nsim = 200;
xref(1,1:Nsim) = 0;
xref(2,1:Nsim) = 0;
xref(3,1:Nsim) = 0.5;
data = loop.simulate(x0, Nsim, 'u.previous', u0, 'x.reference', xref);

%%
figure(1)
plot(data.X(1,:))
legend('x3')
figure(2)
plot(data.X(2,:))
legend('x2')
figure(3)
plot(data.Y)
hold on
plot(xref(3,:))
legend('y')
figure(4)
plot(data.U)
legend('U')