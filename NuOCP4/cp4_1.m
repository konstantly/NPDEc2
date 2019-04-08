%% NuODE 
% CP4
% Konstantinos Brazitikos s1896182
% P 1a 

clear;
close all;

% We make the plots better looking
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesFontSize',20);
fig = gcf; fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;fig.PaperSize = [fig_pos(3) fig_pos(4)];
print('-bestfit','BestFitFigure','-dpdf')

%% delete one %
clf % Delete
format_str = {'interpreter', 'latex','FontSize',30};
KeplerInit;
h1 = 0.001;

totalT = 10;
N = floor(totalT/h1);
[~,Q,~,~]= HamSolver(q0,p0,N,h1,@KeplerForce,'Verlet',Pars);

% Plot for h = 0.001
p0 = plot(0, 0,'bo');
title('Kepler problem with different stepsizes', format_str{:});
xlabel('$x$ position', format_str{:});
ylabel('$y$ position', format_str{:});
axis([-1 2.5 -3.5 1])
hold on
comet(Q(1,:),Q(2,:), 0.3)
legend('Origin','Location', 'SW', format_str{:});
hold on
p1 = plot(Q(1,:),Q(2,:),'r', 'DisplayName', ['$h =$ ', num2str( h1)]);

% Next stepsize  h=0.05
h2 = 0.05; 
KeplerInit;
N = floor(totalT/h2);
[~,Q,~,~]= HamSolver(q0,p0,N,h2,@KeplerForce,'Verlet',Pars);
hold on
p2 = plot(Q(1,:),Q(2,:),'m:', 'DisplayName', ['$h =$ ', num2str( h2)]);

% Next stepsize h = 0.15
h3 = 0.15; 
KeplerInit;
N = floor(totalT/h3);
[~,Q,~,~]= HamSolver(q0,p0,N,h3,@KeplerForce,'Verlet',Pars);
hold on
p3 = plot(Q(1,:),Q(2,:),'c--o', 'DisplayName', ['$h =$ ', num2str( h3)]);

% Next stepsize h = 0.2
h4 = 0.2;
KeplerInit;
N = floor(totalT/h4);
[~,Q,~,~]= HamSolver(q0,p0,N,h4,@KeplerForce,'Verlet',Pars);
hold on
p4 = plot(Q(1,:),Q(2,:),'g--v', 'DisplayName', ['$h =$ ', num2str( h4)]);
%legend([p1 p2 p3 p4 ],{['$h =$ ', num2str( h1)],['$h =$ ' num2str(h2)],...
%    ['$h =$ ' num2str(h3)],['$h =$ ' num2str(h4)]},'Location', 'SW', format_str{:});
hold off
% DELETE
% hold on 
% h = 0.01;
% KeplerInit;
% N = floor(totalT/h);
% [~,Q,~,H1]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
% plot(Q(1,:),Q(2,:),'go')
% 
% KeplerInit;
% h=0.1;
% N = floor(totalT/h);
% [~,Q,~,H1]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
% plot(Q(1,:),Q(2,:),'r-')
% plot(0,0,'*');
% %set(h,'Fontsize',10);
% % xlim([-1,1.5]);
% % ylim([-1.5,1]);
% hold off

%% P 1b
clear;clf;
%%
% We make the plots better looking
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesFontSize',20);
fig = gcf; fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;fig.PaperSize = [fig_pos(3) fig_pos(4)];
%print('-bestfit','BestFitFigure','-dpdf')
format_str = {'interpreter', 'latex','FontSize',30};

KeplerInit;
h=0.001;
totalT = 10;
N = floor(totalT/h);
[~,Q,~,~]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'r-')
hold on
p0=[0.6;0.3];
[~,Q,~,~]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'k-')
p0=[0.6;0.1];
[~,Q,~,~]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'g-')
plot(0,0,'b*')
legend('$\alpha = 0.6$','$\alpha = 0.3$','$\alpha = 0.1$','Fixed point',...
    format_str{:});
xlabel('$x$ position', format_str{:});
ylabel('$y$ position', format_str{:});
title('Kepler problem for different $\alpha$ with $h=0.001$',...
    format_str{:});
%% p2b smaller stepsize h = 0.0001
clear;close all;
set(0,'DefaultLineLineWidth',2); set(0,'DefaultAxesFontSize',20);
KeplerInit;
h=0.0001;
totalT = 10;
N = floor(totalT/h);
[~,Q,~,~]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'r-')
hold on
p0=[0.6;0.3];
[~,Q,~,~]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'k-')
p0=[0.6;0.1];
[~,Q,~,~]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'g-')
plot(0,0,'b*')
h=legend('\alpha = 0.6','\alpha = 0.3','\alpha = 0.1','fixed point')
set(h,'Fontsize',15);
xlabel('x')
ylabel('y')
title('the solution of different \alpha with h=0.0001')
%% p2a b
clear;close all;
set(0,'DefaultLineLineWidth',2); set(0,'DefaultAxesFontSize',20);
KeplerInit;
h=0.001;
totalT = 10;
N = floor(totalT/h);
[~,Q,~,~]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'k^')
hold on
TFCInit;
[~,Q,~,~]= HamSolver(q0,p0,N,h,@TFCForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'r-')
plot(0,0,'b*')
g=legend('KeplerForce','TFCForce','fixed point')
set(g,'Fontsize',15);
xlabel('x')
ylabel('y')
title('the solution of FTC and Kepler')

%% 2c
clear;clf;clc;
axis([-1.5,1.5,-1.5,1.5]);
ax = gca; % current axes
ax.FontSize = 8;
TFCInit;
h = 0.001;
%h=0.0001;
totalT = 10;
N = floor(totalT/h);
[~,Q,~,H_c1]= HamSolver(q0,p0,10000,0.001,@TFCForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'k-')
hold on 
%h=0.00001;
%N = floor(totalT/h);
[~,Q,~,H_c2]= HamSolver(q0,p0,100000,0.0001,@TFCForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'r-')
legend('h= 0.001','h=0.0001')
xlabel('x')
ylabel('y')
title('the solution of different stepsize')
hold off


%%
clear;clf;clc;
ax = gca; % current axes
ax.FontSize = 8;
TFCInit;
h=0.0001;
totalT = 20;
totalT = 10;
N = floor(totalT/h);
[T_verlet,Q,P,H_verlet]= HamSolver(q0,p0,N,h,@TFCForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'k-','LineWidth',2)
hold on
h=0.0004;
N = floor(totalT/h);
[T_rk4,Q,P,H_rk4]= HamSolver(q0,p0,N,h,@TFCForce,'RK4',Pars);
axis([-1.5,1.5,-1.5,1.5]);
plot(Q(1,:),Q(2,:),'r-','LineWidth',2)
legend('h= 0.0001,Verlet','h=0.0004,Rk4')
xlabel('x')
ylabel('y')
title('the solution of different stepsize and methods')
%%
ax = gca; % current axes
ax.FontSize = 8;
subplot(1,3,1)
plot(T_verlet,H_verlet,'k','LineWidth',2)
hold on
plot(T_rk4,H_rk4,'r','LineWidth',2)
axis([0 20 -1.9 -1.5]);
xlabel('time T')
ylabel('energy H')
t = legend('verlet method','RK4 method');
set(t,'Fontsize',6);
subplot(1,3,2)
plot(T_verlet,H_verlet,'k','LineWidth',2)
hold on
plot(T_rk4,H_rk4,'r')
axis([1.35 1.45 -1.875 -1.865]);
xlabel('time T')
ylabel('energy H')
t=legend('verlet method','RK4 method');
set(t,'Fontsize',6);
subplot(1,3,3)
plot(T_verlet,H_verlet,'k','LineWidth',2)
hold on
plot(T_rk4,H_rk4,'r','LineWidth',2)
axis([6 8 -1.875 -1.865]);
xlabel('time T')
ylabel('energy H')
t=legend('verlet method','RK4 method');
set(t,'Fontsize',6);
zoom xon
