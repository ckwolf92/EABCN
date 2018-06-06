figure(1)

subplot(2,3,1)
set(gca,'FontSize',18)
hold on
plot(vTime,c_e_m(1:T_plot+1),'linewidth',2,'linestyle','-','color',[204/255,102/255,0/255])
set(gcf,'color','w')
title('Consumption','interpreter','latex','fontsize',14)
ylabel('$\%$ deviation','interpreter','latex')
xlabel('Horizon','interpreter','latex','FontSize',18)
grid on
hold off

subplot(2,3,2)
set(gca,'FontSize',18)
hold on
plot(vTime,l_e_m(1:T_plot+1),'linewidth',2,'linestyle','-','color',[204/255,102/255,0/255])
set(gcf,'color','w')
title('Labor','interpreter','latex','fontsize',14)
ylabel('$\%$ deviation','interpreter','latex')
xlabel('Horizon','interpreter','latex','FontSize',18)
grid on
hold off

subplot(2,3,3)
set(gca,'FontSize',18)
hold on
plot(vTime,w_e_m(1:T_plot+1),'linewidth',2,'linestyle','-','color',[204/255,102/255,0/255])
set(gcf,'color','w')
title('Wage','interpreter','latex','fontsize',14)
ylabel('$\%$ deviation','interpreter','latex')
xlabel('Horizon','interpreter','latex','FontSize',18)
grid on
hold off

subplot(2,3,4)
set(gca,'FontSize',18)
hold on
plot(vTime,r_n_e_m(1:T_plot+1),'linewidth',2,'linestyle','-','color',[204/255,102/255,0/255])
set(gcf,'color','w')
title('Interest Rate','interpreter','latex','fontsize',14)
ylabel('$\%$ deviation','interpreter','latex')
xlabel('Horizon','interpreter','latex','FontSize',18)
grid on
hold off

subplot(2,3,5)
set(gca,'FontSize',18)
hold on
plot(vTime,y_e_m(1:T_plot+1),'linewidth',2,'linestyle','-','color',[204/255,102/255,0/255])
set(gcf,'color','w')
title('Output','interpreter','latex','fontsize',14)
ylabel('$\%$ deviation','interpreter','latex')
xlabel('Horizon','interpreter','latex','FontSize',18)
grid on
hold off

subplot(2,3,6)
set(gca,'FontSize',18)
hold on
plot(vTime,wedge_EE_e_m(1:T_plot+1),'linewidth',2,'linestyle','-','color',[204/255,102/255,0/255])
set(gcf,'color','w')
title('Wedge','interpreter','latex','fontsize',14)
ylabel('$\%$ deviation','interpreter','latex')
xlabel('Horizon','interpreter','latex','FontSize',18)
grid on
hold off