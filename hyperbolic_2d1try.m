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
N = 100;
T = 100;
dt = 0.1;

% The initial condition
%Pxy as function?
psi_nm1 = exp(-((x-0.25).^2+(y-0.6).^2)/0.08).*chi;

% First "timestep"
% psi_h = psi_nm1 - (dt/2)*result;
psi_h = psi_nm1 - (dt/2)*Jpp(chi, psi_nm1, dx);
psi_n = psi_nm1 - dt.*Jpp(chi, psi_h, dx);

%Do we need this?
psi_np1 = zeros(size(psi_n));
% Remaining "timesteps"
for n = 1:N - 1
    Jn = Jpp(chi, psi_n, dx);
    psi_np1(2:end-1,2:end-1) = psi_nm1(2:end-1,2:end-1)-2*dt*Jn(2:end-1,2:end-1);
    
    psi_nm1 = psi_n;
    psi_n = psi_np1;
    contourf(psi_n)
    drawnow
end

% For two-dimensional plotting you may wish to use contourf or pcolor
