
figure(5)
hold on
for i=1:size(new_Regions,1)
    i
  if i~=12%% && i~=14
        plotregion(-new_Regions{i,1},-new_Regions{i,2}) 
   end
end
plotregion(-A_CR0,-b_CR0)
xlim([-1.5 1.5])
ylim([-1.5 1.5])
% 
% %%
% figure(15)
% x = -1.5:0.01:1.5
% b = new_Regions{12,2};
% A = new_Regions{12,1};
% % i=1;
% % hue = ( b(i) - A(i,1)*x)/A(i,2);
% % %%
% for i = 1:size(b,1)
%     a1x+a2y <= b
%     a2*y <= b - a1*x
%     y <= (b - a1*x)/a2
%     
%     y(i,:) = ( b(i) - A(i,1)*x)/A(i,2); 
%     plot(x,y(i,:))
%     hold on
% end
% 

%%
figure(21)
hold on
for i = 1:5
    A_plot = [-A(i,:); A(6:9,:)]; 
    b_plot = [-b(i,:); b(6:9,:)];
    plotregion(-A_plot,-b_plot)
    xlim([-1.5 1.5])
    ylim([-1.5 1.5])
end
%%
figure(22)
hold on
i=1;
A_plot = [-A(i,:); A(6:9,:)]; 
b_plot = [-b(i,:); b(6:9,:)];
plotregion(-A_plot,-b_plot)
xlim([-1.5 1.5])
ylim([-1.5 1.5])
    %%
A_plot = [A(i,:); A(6:9,:)]; 
b_plot = [b(i,:); b(6:9,:)];
plotregion(-A_plot,-b_plot)
xlim([-1.5 1.5])
ylim([-1.5 1.5])


