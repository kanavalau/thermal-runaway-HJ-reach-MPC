clear
close all

avoid_410 = load('../avoid-MPC/plotting_data_MPC_avoid_410.mat');
stndrd_400 = load('../standard-MPC/plotting_data_MPC_standard_400.mat');
stndrd_410 = load('../standard-MPC/plotting_data_MPC_standard_410.mat');

figure(1)
plot(avoid_410.Sols(1,:),avoid_410.Sols(4,:),'k','LineWidth',1.5)
hold on
plot(stndrd_400.Sols(1,:),stndrd_400.Sols(4,:),'b','LineWidth',1.5)
plot(stndrd_410.Sols(1,:),stndrd_410.Sols(4,:),'r','LineWidth',1.5)
axis([0 5000 395 450])
ax = gca;
ax.FontSize = 14;
ax.TickLabelInterpreter = 'latex';
xlabel('Time, $t$/s','interpreter','latex','FontSize',14)
ylabel('Reactor temperature $T_\mathrm{r}$/K','interpreter','latex','FontSize',14)
legend('Avoid-MPC $T_\mathrm{sp}=410$ K','Standard-MPC $T_\mathrm{sp}=400$ K','Standard-MPC $T_\mathrm{sp}=410$ K',...
    'interpreter','latex','FontSize',14)

figure(2)
plot(avoid_410.Sols(1,:),avoid_410.Sols(3,:),'k','LineWidth',1.5)
hold on
plot(stndrd_400.Sols(1,:),stndrd_400.Sols(3,:),'b','LineWidth',1.5)
plot(stndrd_410.Sols(1,:),stndrd_410.Sols(3,:),'r','LineWidth',1.5)
axis([0 5000 0 8])
ax = gca;
ax.FontSize = 14;
ax.TickLabelInterpreter = 'latex';
xlabel('Time, $t$/s','interpreter','latex','FontSize',14)
ylabel('Product concentration $c_\mathrm{B}$/kmol m$^{-3}$','interpreter','latex','FontSize',14)
legend('Avoid-MPC $T_\mathrm{sp}=410$ K','Standard-MPC $T_\mathrm{sp}=400$ K','Standard-MPC $T_\mathrm{sp}=410$ K',...
    'interpreter','latex','FontSize',14)

figure(3)
plot(avoid_410.Sols(1,:),avoid_410.alpha,'k','LineWidth',1.5)
hold on
plot(stndrd_400.Sols(1,:),stndrd_400.alpha,'b','LineWidth',1.5)
plot(stndrd_410.Sols(1,1:2105),stndrd_410.alpha(1:2105),'r','LineWidth',1.5)
ax = gca;
ax.FontSize = 24;
ax.TickLabelInterpreter = 'latex';
xlabel('Time, $t$/s','interpreter','latex','FontSize',24)
ylabel('Level-set function $\alpha^{*}$','interpreter','latex','FontSize',24)
%legend('Avoid-MPC $T_\mathrm{sp}=410$ K','Standard-MPC $T_\mathrm{sp}=400$ K','Standard-MPC $T_\mathrm{sp}=410$ K',...
%    'interpreter','latex','FontSize',20,'Location','SouthEast')

figure(4)
plot(avoid_410.Sols(1,:),avoid_410.alpha,'k','LineWidth',1.5)
hold on
plot(stndrd_400.Sols(1,:),stndrd_400.alpha,'b','LineWidth',1.5)
plot(stndrd_410.Sols(1,1:2105),stndrd_410.alpha(1:2105),'r','LineWidth',1.5)
axis([0 5000 -40 60])
ax = gca;
ax.FontSize = 24;
ax.TickLabelInterpreter = 'latex';
xlabel('Time, $t$/s','interpreter','latex','FontSize',24)
ylabel('Level-set function $\alpha^{*}$','interpreter','latex','FontSize',24)
legend('Avoid-MPC $T_\mathrm{sp}=410$ K','Standard-MPC $T_\mathrm{sp}=400$ K','Standard-MPC $T_\mathrm{sp}=410$ K',...
    'interpreter','latex','FontSize',24,'Location','SouthEast')