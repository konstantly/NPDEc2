clear all;
%close all;

% Define a "spatial" grid
M = 100;
D = 1;
x = linspace(0, D, M + 1);
y = linspace(0, D, M + 1);
[X, Y] = meshgrid(x, y);
X = X';
Y = Y';
% X(m, p) = x_{m - 1} for all p
% Y(m, p) = y_{p - 1} for all m
dx = D / M;

chi = sin(pi*X).*sin(pi*Y);

% Define information about the step sizes in the "time" dimension
% Not sure about those values
T = 10;
c = 5;
dt = dx/c; %Is this ok?
N = floor(T * (1.0 + 1.0e-10) / dt)

% The initial condition
%Pxy as function?
psi_nm1 = exp(-((x-0.25).^2+(y-0.6).^2)/0.08).*chi;

% First "timestep"
% psi_h = psi_nm1 - (dt/2)*result;

psi_nm11 = psi_nm1;
psi_nm12 = psi_nm1;
psi_nm13 = psi_nm1;

% We use the variant of Jnmp that we want each time
% by commenting in and out the corresponding lines here
% and in the calculation of the next step.
% J++
psi_h1 = psi_nm11 - (dt/2)*Jpp(chi, psi_nm1, dx);
psi_n1 = psi_nm11 - dt.*Jpp(chi, psi_h1, dx);

% J+x
psi_h2 = psi_nm12 - (dt/2)*Jpx(chi, psi_nm1, dx);
psi_n2 = psi_nm12 - dt.*Jpx(chi, psi_h2, dx);

% Jx+
psi_h3 = psi_nm13 - (dt/2)*Jxp(chi, psi_nm1, dx);
psi_n3 = psi_nm13 - dt.*Jxp(chi, psi_h3, dx);

%Do we need this?
psi_np11 = zeros(size(psi_n1));
psi_np12 = zeros(size(psi_n2));
psi_np13 = zeros(size(psi_n3));
t = zeros(1,N-1);
sumPsi1 = zeros(1,N-1);
sumPsi2 = zeros(1,N-1);
sumPsi3 = zeros(1,N-1);
sumSqPsi1 = zeros(1,N-1);
sumSqPsi2 = zeros(1,N-1);
sumSqPsi3 = zeros(1,N-1);


% We make the plots better looking
% set(0, 'DefaultLineLineWidth', 1.2);
% set(0, 'DefaultAxesFontSize',13);
% fig = gcf; fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print('-bestfit','BestFitFigure','-dpdf')
% format_str = {'interpreter', 'latex','FontSize',25};

format_str = {'Interpreter', 'latex','FontSize',30};
set(0, 'DefaultAxesFontSize',20);
set(0,'DefaultLineLineWidth', 2);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
fig = gcf; fig.PaperPositionMode = 'auto'; fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];


% Remaining "timesteps"
for n = 1:N - 1
    
    % We select our Jnmp 
    Jn1 = Jpp(chi, psi_n1, dx);
    Jn2 = Jpx(chi, psi_n2, dx);
    Jn3 = Jxp(chi, psi_n3, dx);
    
    % Evaluate the value at Psi n+1
    psi_np11(2:end-1,2:end-1) = psi_nm11(2:end-1,2:end-1)-2*dt*Jn1(2:end-1,2:end-1);
    psi_np12(2:end-1,2:end-1) = psi_nm12(2:end-1,2:end-1)-2*dt*Jn2(2:end-1,2:end-1);
    psi_np13(2:end-1,2:end-1) = psi_nm13(2:end-1,2:end-1)-2*dt*Jn3(2:end-1,2:end-1);
    
    % reset the step in order to perform the next iteration
    psi_nm11 = psi_n1;
    psi_nm12 = psi_n2;
    psi_nm13 = psi_n3;

    psi_n1 = psi_np11;
    psi_n2 = psi_np12;
    psi_n3 = psi_np13;

    t(n) = n*dt;
    sumPsi1(n) = dx^2*sum(sum(psi_n1));
    sumPsi2(n) = dx^2*sum(sum(psi_n2));
    sumPsi3(n) = dx^2*sum(sum(psi_n3));
    sumSqPsi1(n) = dx^2*sum(sum(psi_n1.*psi_n1));
    sumSqPsi2(n) = dx^2*sum(sum(psi_n2.*psi_n2));
    sumSqPsi3(n) = dx^2*sum(sum(psi_n3.*psi_n3));

%     figure(1)
%     pcolor(psi_n1)
%     contourf(psi_n1)
%     title('Density of $\Psi$ using J++', format_str{:});
%     xlabel('$x$ grid', format_str{:});
%     ylabel('$y$ grid', format_str{:});
%     drawnow()

%     figure(2)
%     pcolor(psi_n2)
%     contourf(psi_n2)
%     title('Density of $\Psi$ using J+x', format_str{:});
%     xlabel('$x$ grid', format_str{:});
%     ylabel('$y$ grid', format_str{:});
%     drawnow()

%     figure(3)
%     pcolor(psi_n3)
%     contourf(psi_n3)
%     title('Density of $\Psi$ using Jx+', format_str{:});
%     xlabel('$x$ grid', format_str{:});
%     ylabel('$y$ grid', format_str{:});
%     drawnow()

    if mod(n, 500) == 0
        disp(n);
    end
end
%%
% Make them with latex
figure(2)
%subplot(1,2,1)
plot(t, sumPsi1)
hold on
plot(t, sumPsi2)
hold on
plot(t, sumPsi3)
hold off
title('Sum of $\psi$ against time.', format_str{:});
xlabel('Time $t$', format_str{:});
ylabel('Sum', format_str{:});
legend('$J^{++}$', '$J^{+\times}$', '$J^{\times+}$',...
    'Location','NW')

%% subplot(1,2,2)
plot(t, sumSqPsi1)
hold on
plot(t, sumSqPsi2)
hold on
plot(t, sumSqPsi3)
hold off
title('Sum of Squares of $\psi$ against time.', format_str{:});
xlabel('Time $t$', format_str{:});
ylabel('Sum', format_str{:});
legend('$J^{++}$', '$J^{+\times}$', '$J^{\times+}$',...
    'Location','NW')
%%
%title('Error of $\Psi$', format_str{:});
%legend('sum of \Psi', 'sum of Squares of \Psi')
%xlabel('$n\Delta t$', format_str{:});
%ylabel('sum of $\Psi$ times $\Delta x$', format_str{:});
%suptitle('I am a super title')
%ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],...
%    'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%text(0.5, 0.98,'Title', format_str{:});
%sgtitle('asd', format_str{:});

% For two-dimensional plotting you may wish to use contourf or pcolor
