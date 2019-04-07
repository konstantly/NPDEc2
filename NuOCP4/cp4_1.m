clear;close all;
set(0,'DefaultLineLineWidth',2); set(0,'DefaultAxesFontSize',20);
KeplerInit;
h = 0.001;
%h=0.1;

N=10/h;
[T,Q,P,H1]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'k^')
hold on 
h = 0.01;
KeplerInit;
N=10/h;
[T,Q,P,H1]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'go')

KeplerInit;
h=0.1;
N=10/h;
[T,Q,P,H1]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'r-')
title('the solution of Kepler with different stepsize')
plot(0,0,'*');
h=legend('h = 0.001','h=0.01','h=0.1','fixedpoint');
set(h,'Fontsize',10);
xlim([-1,1.5]);
ylim([-1.5,1]);
xlabel('x')
ylabel('y')
hold off
%% p1b
clear;close all;
set(0,'DefaultLineLineWidth',2); set(0,'DefaultAxesFontSize',20);
KeplerInit;
h=0.001;
N=10/h;
[T,Q,P,H4]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'r-')
hold on
p0=[0.6;0.3];
[T,Q,P,H]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'k-')
p0=[0.6;0.1];
[T,Q,P,H]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'g-')
plot(0,0,'b*')
legend('\alpha = 0.6','\alpha = 0.3','\alpha = 0.1','fixed point')
xlabel('x')
ylabel('y')
title('the solution of different \alpha with h=0.001')
%% p2b smaller stepsize h = 0.0001
clear;close all;
set(0,'DefaultLineLineWidth',2); set(0,'DefaultAxesFontSize',20);
KeplerInit;
h=0.0001;
N=10/h;
[T,Q,P,H4]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'r-')
hold on
p0=[0.6;0.3];
[T,Q,P,H]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'k-')
p0=[0.6;0.1];
[T,Q,P,H]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
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
N=10/h;
[T,Q,P,H]= HamSolver(q0,p0,N,h,@KeplerForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'k^')
hold on
TFCInit;
[T,Q,P,H]= HamSolver(q0,p0,N,h,@TFCForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'r-')
plot(0,0,'b*')
g=legend('KeplerForce','TFCForce','fixed point')
set(g,'Fontsize',15);
xlabel('x')
ylabel('y')
title('the solution of FTC and Kepler')

%% 2c
clear;close all;clc;
axis([-1.5,1.5,-1.5,1.5]);
ax = gca; % current axes
ax.FontSize = 8;
TFCInit;
h = 0.001;
%h=0.0001;
N=10/h;
[T,Q,P,H_c1]= HamSolver(q0,p0,10000,0.001,@TFCForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'k-')
hold on 
%h=0.00001;
%N=10/h;
[T,Q,P,H_c2]= HamSolver(q0,p0,100000,0.0001,@TFCForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'r-')
legend('h= 0.001','h=0.0001')
xlabel('x')
ylabel('y')
title('the solution of different stepsize')
hold off


%%
clear;close all;clc;
ax = gca; % current axes
ax.FontSize = 8;
TFCInit;
h=0.0001;
N=20/h;
[T_verlet,Q,P,H_verlet]= HamSolver(q0,p0,N,h,@TFCForce,'Verlet',Pars);
plot(Q(1,:),Q(2,:),'k-','LineWidth',2)
hold on
h=0.0004;
N=20/h;
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
