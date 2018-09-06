% %     clear f
% %     clear x
% %     clear LMI
% %     %fk = sdpvar(1,1,'full'); 
% %     xk = sdpvar(2,1,'full'); 
% %     n = 1;
% %     A = A_CR0
% %     b = b_CR0
% %     m = size(A,1);
% %     LMI = [];
% %     %LMI = [LMI, -G*z + W + S*x0 >= 0 ];
% %     LMI = [LMI, A(2:7,:)*xk <=b(2:7)];
% %     LMI = [LMI, A(1,:)*xk <= (b(1)+1)]
% %     objetivo = A(1,:)*xk;
% %     options = sdpsettings;
% %     options.solver = 'sedumi';
% % 
% %     diagnostics = optimize(LMI,-objetivo,options)
% %     
% %     double(objetivo)
% %     b(1)
% %     if (double(objetivo)) <= (b(1))
% %         disp('eh redundante')
% %     end
% %     
%     
% %
% figure(10)
% hold on
% for i=1:size(CRest_3,1)
%     i
%   %if i~=10 && i~=14 && i~=15
%         plotregion(-CRest_3{i,1},-CRest_3{i,2}) 
%   %end
%    xlim([-1.5 1.5])
% ylim([-1.5 1.5])
% end
% plotregion(-A_CR0,-b_CR0)
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])

%%
figure(5)
hold on
for i=1:size(new_Regions,1)
    i
  %if i~=10 && i~=15 %&& i~=15
        plotregion(-new_Regions{i,1},-new_Regions{i,2}) 
  %end
   xlim([-1.5 1.5])
ylim([-1.5 1.5])
end
plotregion(-A_CR0,-b_CR0)
xlim([-1.5 1.5])
ylim([-1.5 1.5])
% 
%%
% figure(17)
% x = -1.5:0.01:1.5
% b = new_Regions{9,2};
% A = new_Regions{9,1};
% % i=1;
% % hue = ( b(i) - A(i,1)*x)/A(i,2);
% % %%
% for i = 1:size(b,1)
%     %a1x+a2y <= b
%     %a2*y <= b - a1*x
%     %y <= (b - a1*x)/a2
%     
%     y(i,:) = ( b(i) - A(i,1)*x)/A(i,2); 
%     plot(x,y(i,:))
%     hold on
% end


% %%
% figure(21)
% hold on
% for i = 1:5
%     A_plot = [-A(i,:); A(6:9,:)]; 
%     b_plot = [-b(i,:); b(6:9,:)];
%     plotregion(-A_plot,-b_plot)
%     xlim([-1.5 1.5])
%     ylim([-1.5 1.5])
% end
% %%
% figure(22)
% hold on
% i=1;
% A_plot = [-A(i,:); A(6:9,:)]; 
% b_plot = [-b(i,:); b(6:9,:)];
% plotregion(-A_plot,-b_plot)
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])
%     %%
% A_plot = [A(i,:); A(6:9,:)]; 
% b_plot = [b(i,:); b(6:9,:)];
% plotregion(-A_plot,-b_plot)
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])
% 
% 
