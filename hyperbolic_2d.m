clear all;
close all;

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
c = 2;
dt = dx/c; %Is this ok?
N = floor(T * (1.0 + 1.0e-10) / dt)

% The initial condition
%Pxy as function?
psi_nm1 = exp(-((x-0.25).^2+(y-0.6).^2)/0.08).*chi;

% First "timestep"
% psi_h = psi_nm1 - (dt/2)*result;

% We use the variant of Jnmp that we want each time
% by commenting in and out the corresponding lines here
% and in the calculation of the next step.
% J++
% psi_h = psi_nm1 - (dt/2)*Jpp(chi, psi_nm1, dx);
% psi_n = psi_nm1 - dt.*Jpp(chi, psi_h, dx);

% J+x
psi_h = psi_nm1 - (dt/2)*Jpx(chi, psi_nm1, dx);
psi_n = psi_nm1 - dt.*Jpx(chi, psi_h, dx);

% Jx+
% psi_h = psi_nm1 - (dt/2)*Jxp(chi, psi_nm1, dx);
% psi_n = psi_nm1 - dt.*Jxp(chi, psi_h, dx);

%Do we need this?
psi_np1 = zeros(size(psi_n));

% We make the plots better looking
set(0, 'DefaultLineLineWidth', 1.2);
set(0, 'DefaultAxesFontSize',13);
fig = gcf; fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;fig.PaperSize = [fig_pos(3) fig_pos(4)];
print('-bestfit','BestFitFigure','-dpdf')
format_str = {'interpreter', 'latex','FontSize',25};

% Remaining "timesteps"
for n = 1:N - 1
    
    % We select our Jnmp 
%     Jn = Jpp(chi, psi_n, dx);
    Jn = Jpx(chi, psi_n, dx);
%     Jn = Jxp(chi, psi_n, dx);
    
    % Evaluate the value at Psi n+1
    psi_np1(2:end-1,2:end-1) = psi_nm1(2:end-1,2:end-1)-2*dt*Jn(2:end-1,2:end-1);
    
    % reset the step in order to perform the next iteration
    psi_nm1 = psi_n;
    psi_n = psi_np1;
    t(n) = n*dt;
    sumPsi(n) = dx^2*sum(sum(psi_n));
    sumSqPsi(n) = dx^2*sum(sum(psi_n.*psi_n));
    
%    pcolor(psi_n)
    contourf(psi_n)
    title('Density of $\Psi$', format_str{:});
    xlabel('$x$ grid', format_str{:});
    ylabel('$y$ grid', format_str{:});
    drawnow()
end
%%
% Make them with latex
figure(2)
plot(t, sumPsi)
hold on
plot(t, sumSqPsi)
title('Error of $\Psi$', format_str{:});
legend('sum of \Psi', 'sum of Squares of \Psi')
xlabel('$n\Delta t$', format_str{:});
ylabel('sum of $\Psi$ times $\Delta x$', format_str{:});
hold off
% For two-dimensional plotting you may wish to use contourf or pcolor
