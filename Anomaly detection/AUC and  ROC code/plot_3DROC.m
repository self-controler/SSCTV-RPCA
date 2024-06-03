function plot_3DROC(det_map,GT,detec_label,mode_eq)
% det_map: the detection result N*k, N is the number of pixels in the detection result, 
% GT: Ground truth
% detec_label: the name of detector for legend
% mode_eq: if mode_eq==1. the equation (7) in the paper is used;the equation (9)is used;

num_map = size(det_map,2);
for i = 1:num_map
    det_map(:,i) = (det_map(:,i) - min(det_map(:,i))) /(max(det_map(:,i))-min(det_map(:,i)));
end

%PD and PF based on uniform step and sample value
for k = 1:num_map
    tau1(:,k) = [0:0.01:1]';
end
 tau1 = sort(tau1,'descend');

for k = 1:num_map
    for i = 1: length(tau1)
        map =det_map(:,k);
        if mode_eq==1
           map(det_map(:,k)>=tau1(i,k))=1;
           map(det_map(:,k)<tau1(i,k))=0;
        else
           map(det_map(:,k)>tau1(i,k))=1;
           map(det_map(:,k)<=tau1(i,k))=0;
        end
        [PD1(i,k),PF1(i,k)] = cal_pdpf(map,GT);
    end
end
map = [];

tau2 = sort(det_map,'descend');

 
 for k = 1:num_map
    for i = 1: length(tau2)
        map =det_map(:,k);
        if mode_eq==1
           map(det_map(:,k)>=tau2(i,k))=1;
           map(det_map(:,k)<tau2(i,k))=0;
        else
           map(det_map(:,k)>tau2(i,k))=1;
           map(det_map(:,k)<=tau2(i,k))=0;
        end
        [PD2(i,k),PF2(i,k)] = cal_pdpf(map,GT);
    end
 end
 
% 
 a11 = min(PD1(1,:));
 a10 = max(max(PD1)); 
 b11=min(PF1(1,:));
 b10 = max(max(PF1));
 
 a21 = min(PD2(1,:));
 a20 = max(max(PD2)); 
 b21=min(PF2(1,:));
 b20 = max(max(PF2)); 
PD1nor = (PD1-a11)/(a10-a11);
PF1nor = (PF1-b11)/(b10-b11);
PD2nor = (PD2-a21)/(a20-a21);
PF2nor = (PF2-b21)/(b20-b21);
 
% % show ROC (PF, PD)
% figure,
% %plot(PF1,PD1,'LineWidth',2)
% 
% semilogx((PF1(:,1))',(PD1(:,1))','Color',[0.5,0.16,0.16],'LineWidth',1.5);
% hold on;
% semilogx((PF1(:,2))',(PD1(:,2))','-k','LineWidth',1.5);
% semilogx((PF1(:,3))',(PD1(:,3))','-y','LineWidth',1.5);
% semilogx((PF1(:,4))',(PD1(:,4))','-m','LineWidth',1.5);
% semilogx((PF1(:,5))',(PD1(:,5))','-c','LineWidth',1.5);
% semilogx((PF1(:,6))',(PD1(:,6))','color',[1,0.5,0],'LineWidth',1.5);
% semilogx((PF1(:,7))',(PD1(:,7))','-g','LineWidth',1.5);
% semilogx((PF1(:,8))',(PD1(:,8))','-b','LineWidth',1.5);
% semilogx((PF1(:,9))',(PD1(:,9))','-r','LineWidth',1.5);
% 
% 
% set(gca,'XTick',(0:0.2:1),'fontsize',16)
% set(gca,'YTick',(0:0.2:1),'fontsize',16)
% %set(gca, 'XDir','reverse')
% xlabel('P_F','Fontname', 'Times New Roman','fontsize',18) ; ylabel('P_D','Fontname', 'Times New Roman','fontsize',18) 
% grid on

% show ROC (PF, PD)
figure,
%plot(PF1,PD1,'LineWidth',2)

plot((PF1(:,1))',(PD1(:,1))','--','Color',[0.5,0.16,0.16],'LineWidth',1.5);
hold on;
plot((PF1(:,2))',(PD1(:,2))','--k','LineWidth',1.5);
plot((PF1(:,3))',(PD1(:,3))','--y','LineWidth',1.5);
plot((PF1(:,4))',(PD1(:,4))','--m','LineWidth',1.5);
plot((PF1(:,5))',(PD1(:,5))','--c','LineWidth',1.5);
plot((PF1(:,6))',(PD1(:,6))','--','color',[1,0.5,0],'LineWidth',1.5);
plot((PF1(:,7))',(PD1(:,7))','-g','LineWidth',1.5);
plot((PF1(:,8))',(PD1(:,8))','-b','LineWidth',1.5);
plot((PF1(:,9))',(PD1(:,9))','-r','LineWidth',1.5);

axis([1e-4 1 0 1])
set(gca,'xscale','log')
%set(gca,'XTick',(0:0.2:1),'fontsize',16)
set(gca,'YTick',(0:0.2:1),'fontsize',16)
%set(gca, 'XDir','reverse')

xlabel('P_F','Fontname', 'Times New Roman','fontsize',18) ; ylabel('P_D','Fontname', 'Times New Roman','fontsize',18) 
grid on


for i = 1:num_map
    name1(i) =strcat(detec_label(i)); 
    name2(i) = strcat(detec_label(i));
end
legend([name1,name2],'Location','southeast')
%legend boxoff


% show ROC (PD, Tau) based on uniform step and sample value

figure,

plot((tau1(:,1))',(PD1(:,1))','--','Color',[0.5,0.16,0.16],'LineWidth',1.5);
hold on;
plot((tau1(:,2))',(PD1(:,2))','--k','LineWidth',1.5);
plot((tau1(:,3))',(PD1(:,3))','--y','LineWidth',1.5);
plot((tau1(:,4))',(PD1(:,4))','--m','LineWidth',1.5);
plot((tau1(:,5))',(PD1(:,5))','--c','LineWidth',1.5);
plot((tau1(:,6))',(PD1(:,6))','--','color',[1,0.5,0],'LineWidth',1.5);
plot((tau1(:,7))',(PD1(:,7))','-g','LineWidth',1.5);
plot((tau1(:,8))',(PD1(:,8))','-b','LineWidth',1.5);
plot((tau1(:,9))',(PD1(:,9))','-r','LineWidth',1.5);
axis([0,1,0,1])

%axis([1e-4 1 0 1])
set(gca,'xscale','log')
%set(gca,'XTick',(0:0.2:1),'fontsize',16)
set(gca,'YTick',(0:0.2:1),'fontsize',16)
xlabel('\tau','Fontname', 'Times New Roman','fontsize',18) ; 

ylabel('P_D','Fontname', 'Times New Roman','fontsize',18) 
grid on

%set(gca, 'XDir','reverse')

hold off
legend([name1,name2],'Location','southeast')
%legend boxoff


% show ROC (PF, Tau) based on uniform step and sample value
figure,
%plot(tau2,PF2,'LineWidth',2,'linestyle','--')

plot((tau1(:,1))',(PF1(:,1))','--','Color',[0.5,0.16,0.16],'LineWidth',1.5);
hold on;
plot((tau1(:,2))',(PF1(:,2))','--k','LineWidth',1.5);
plot((tau1(:,3))',(PF1(:,3))','--y','LineWidth',1.5);
plot((tau1(:,4))',(PF1(:,4))','--m','LineWidth',1.5);
plot((tau1(:,5))',(PF1(:,5))','--c','LineWidth',1.5);
plot((tau1(:,6))',(PF1(:,6))','--','color',[1,0.5,0],'LineWidth',1.5);
plot((tau1(:,7))',(PF1(:,7))','-g','LineWidth',1.5);
plot((tau1(:,8))',(PF1(:,8))','-b','LineWidth',1.5);
plot((tau1(:,9))',(PF1(:,9))','-r','LineWidth',1.5);
axis([0,1,0,1])

%axis([1e-3 1 0 1])
set(gca,'xscale','log')
%set(gca,'XTick',(0:0.2:1),'fontsize',16)
set(gca,'YTick',(0:0.2:1),'fontsize',16)
xlabel('\tau','Fontname', 'Times New Roman','fontsize',18) ; ylabel('P_F','Fontname', 'Times New Roman','fontsize',18) 
grid on
 %semilogx(FPR(k), TPR(k), 'bo','linewidth',2.5,'MarkerSize',4);

hold off
legend([name1,name2],'Location','southeast');
%legend boxoff



% 3D ROC

figure,%plot3(PF1,tau1,PD1,'LineWidth',2)

plot3(PF1(:,1),tau1(:,1),PD1(:,1),'--','Color',[0.5,0.16,0.16],'LineWidth',1.5);
hold on;
plot3(PF1(:,2),tau1(:,2),PD1(:,2),'--k','LineWidth',1.5);
plot3(PF1(:,3),tau1(:,3),PD1(:,3),'--y','LineWidth',1.5);
plot3(PF1(:,4),tau1(:,4),PD1(:,4),'--m','LineWidth',1.5);
plot3(PF1(:,5),tau1(:,5),PD1(:,5),'--c','LineWidth',1.5);
plot3(PF1(:,6),tau1(:,6),PD1(:,6),'--','color',[1,0.5,0],'LineWidth',1.5);
plot3(PF1(:,7),tau1(:,7),PD1(:,7),'-g','LineWidth',1.5);
plot3(PF1(:,8),tau1(:,8),PD1(:,8),'-b','LineWidth',1.5);
plot3(PF1(:,9),tau1(:,9),PD1(:,9),'-r','LineWidth',1.5);


axis([0, 1, 0, 1, 0, 1])
set(gca,'XTick',(0:0.2:1),'fontsize',16)
set(gca,'YTick',(0:0.2:1),'fontsize',16)
set(gca,'ZTick',(0:0.2:1),'fontsize',16)
xlabel('P_F','Fontname', 'Times New Roman','fontsize',18); 
ylabel('\tau','Fontname', 'Times New Roman','fontsize',18); 
zlabel('P_D','Fontname','Times New Roman','fontsize',18); 
grid on

ax = gca;
ax.BoxStyle = 'full';
box on
legend([name1,name2],'Location','southeast')
%legend boxoff


