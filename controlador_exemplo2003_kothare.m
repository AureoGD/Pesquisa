close all

figure(1)
hold on
for i=1:S
    [x(:,S+1-i) y(:,S+1-i)] = plot_elipse(inQ_table{i},1000);
    F1(S+1-i) = F_table{1,i}(1);
    F2(S+1-i) = F_table{1,i}(2);
    X1_dec(i) = X1(1,S+1-i);
end


figure(2)
semilogx(X1(1,:),F1)

figure(3)
semilogx(X1(1,:),F2)

%%
for i = 2:S
    A = double(Q_table{i-1}) - double(Q_table{i});
    [lixo, p] = chol(A,'lower');
    if p ~=0
        i
    end
end

%%
figure(4)
hold on
for i=1:S
    if i == 1 || i == 12 || i == 13
        plot(y(:,i),-x(:,i))
    else
        plot(x(:,i),y(:,i))
    end
 end