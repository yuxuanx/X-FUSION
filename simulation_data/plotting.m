load('10_98')
% errorbar(x,y,yneg,ypos,xneg,xpos,'o','MarkerSize',5,...
%     'MarkerEdgeColor','red','MarkerFaceColor','red','Linewidth',1);
% figure('units','normalized','position',[.2 .2 .4 .4])
figure(1)
subplot(2,2,1)
for i = 1:9
    errorbar(x(i),y1(i),yneg1(i),ypos1(i),xneg(i),xpos(i),'o','MarkerSize',5,'Linewidth',1);
    hold on
end
set(gca,'xscale','log')
xlabel('Time (seconds)')
ylabel('MGOSPA (Total)')
legend('\delta-GLMB','\delta-GLMB(Joint)','LMB','PMBM','PMBM(Recycle)',...
    'PMB(LBP)','PMB(LBP+Recycle)','PMB(Murty)','PMB(Murty+Recycle)','Location','best');
subplot(2,2,2)
for i = 1:9
    errorbar(x(i),y2(i),yneg2(i),ypos2(i),xneg(i),xpos(i),'o','MarkerSize',5,'Linewidth',1);
    hold on
end
set(gca,'xscale','log')
xlabel('Time (seconds)')
ylabel('MGOSPA (Loc)')
% legend('\delta-GLMB','\delta-GLMB(Joint)','LMB','PMBM','PMBM(Recycle)',...
%     'PMB(LBP)','PMB(LBP+Recycle)','PMB(Murty)','PMB(Murty+Recycle)','Location','best');
subplot(2,2,3)
for i = 1:9
    errorbar(x(i),y3(i),yneg3(i),ypos3(i),xneg(i),xpos(i),'o','MarkerSize',5,'Linewidth',1);
    hold on
end
set(gca,'xscale','log')
xlabel('Time (seconds)')
ylabel('MGOSPA (Missed)')
% legend('\delta-GLMB','\delta-GLMB(Joint)','LMB','PMBM','PMBM(Recycle)',...
%     'PMB(LBP)','PMB(LBP+Recycle)','PMB(Murty)','PMB(Murty+Recycle)','Location','best');
subplot(2,2,4)
for i = 1:9
    errorbar(x(i),y4(i),yneg4(i),ypos4(i),xneg(i),xpos(i),'o','MarkerSize',5,'Linewidth',1);
    hold on
end
set(gca,'xscale','log')
xlabel('Time (seconds)')
ylabel('MGOSPA (False)')
% legend('\delta-GLMB','\delta-GLMB(Joint)','LMB','PMBM','PMBM(Recycle)',...
%     'PMB(LBP)','PMB(LBP+Recycle)','PMB(Murty)','PMB(Murty+Recycle)','Location','best');

