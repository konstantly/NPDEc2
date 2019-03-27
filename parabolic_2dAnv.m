clear all;
close all;

% Define a "spatial" grid
M = 40:20:200;
dx_list = zeros(1, length(M));
err_norm1_list = zeros(1, length(M));
err_norm2_list = zeros(1, length(M));
er_norm1 = 0;
er_norm2 = 0;
for k=1:length(M)
    D = 1;
    x = linspace(0, D, M(k) + 1);
    y = linspace(0, D, M(k) + 1);
    [X, Y] = meshgrid(x, y);
    X = X';
    Y = Y';
    % X(m, p) = x_{m - 1} for all p
    % Y(m, p) = y_{p - 1} for all m
    dx = D / M(k);
    dx_list(k) = dx;
    

    % Define information about the step sizes in the "time" dimension
    dt = dx^2/8;
    T = 0.1;
    N = floor(T * (1.0 + 1.0e-10) / dt);

    alpha = 1;

    % The initial condition
    psi_n = 4*sin(4*pi*X).*sin(5*pi*Y)+sin(pi*X).*sin(pi*Y);

    psi_np1 = zeros(size(psi_n));

    er_norm1 = 0;
    er_norm2 = 0;
    %error_norm = 

    for n = 0:N - 1
    psi_np1(2:end-1,2:end-1) = psi_n(2:end-1,2:end-1) - dt*psi_n(2:end-1,2:end-1)...
        + (dt/(dx^2))*(1-alpha)*(psi_n(2:end-1,3:end)+psi_n(1:end-2,2:end-1)...
        - 4*psi_n(2:end-1,2:end-1) + psi_n(3:end,2:end-1)+ psi_n(2:end-1,1:end-2))...
        + (dt/(2*dx^2))*alpha*(psi_n(1:end-2,3:end)+psi_n(3:end,3:end)-4*psi_n(2:end-1,2:end-1)...
        + psi_n(1:end-2,1:end-2) + psi_n(3:end,1:end-2));

  
    % Exact solution after (n+1) steps
    Psi_np1 = 4*exp(-dt*(n+1)-(16+25)*pi^2*(n+1)*dt).*sin(4*pi*X).*sin(5*pi*Y) ...
        + exp(-dt*(n+1)-(1+1)*(n+1)*pi^2*dt).*sin(pi*X).*sin(pi*Y);
  
    err = (psi_np1 - Psi_np1);
    
    er_norm1 = max(er_norm1 ,sqrt(sum(sum(err.^2))));
    er_norm2 = max(er_norm2, sqrt(dx*dx.*sum(sum(err.^2))));

    % Plot the solution on the (n + 1)th timestep
%     figure(1)
%     plot(y, psi_np1, 'k-')
%     xlabel('$x$', 'FontSize', 16, 'Interpreter', 'latex')
%     ylabel('Numerical solution', 'FontSize', 16)
%     xlim([0.0, 1.0])
%     text(0.55, 2, ['$t = ', num2str((n + 1) * dt, '%.4f'), ', n = ', num2str(n + 1), '$'], 'FontSize', 16, 'Interpreter', 'latex')
%     drawnow()

    psi_n = psi_np1;
    
    figure(2)
    contourf(psi_n)
%     figure(3)
%     plot3(X, Y, psi_np1, 'k-')
%     figure(4)
%     contour3(X,Y,psi_np1)
%     figure(5)
%     surf(X,Y,psi_np1)
    drawnow()
    end
    err_norm1_list(k) = abs(max(er_norm1));
    err_norm2_list(k) = abs(max(er_norm2));

% For two-dimensional plotting you may wish to use contourf or pcolor
end
figure(6)
loglog(dx_list, err_norm1_list)
hold on
loglog(dx_list, err_norm2_list)
hold off
xlabel({'$\Delta x$'}, 'Interpreter', 'latex')
ylabel('Error')
title('Error Vs \Delta x')  
legend({'Error of 2 norm', 'Error of 2 \Delta x norm'})
slope =(log(err_norm1_list(2))-log(err_norm1_list(1)))/(log(dx_list(2))-log(dx_list(1)));
fprintf('The slope for the 2 norm comes out to be: %d \n', slope);
slope =(log(err_norm2_list(2))-log(err_norm2_list(1)))/(log(dx_list(2))-log(dx_list(1)));
fprintf('The slope for the 2 dx norm comes out to be: %d \n', slope);

%%

clear all
close all

% Define a "spatial" grid
M = 40:20:200;
dx_list = zeros(1, length(M));
err_norm1_list = zeros(1, length(M));
err_norm2_list = zeros(1, length(M));
er_norm1 = 0;
er_norm2 = 0;
ratio = linspace(0.01,1,length(M));
for k=1:length(M)
    D = 1;
    x = linspace(0, D, M(k) + 1);
    y = linspace(0, D, M(k) + 1);
    [X, Y] = meshgrid(x, y);
    X = X';
    Y = Y';
    % X(m, p) = x_{m - 1} for all p
    % Y(m, p) = y_{p - 1} for all m
    dx = D / M(k);
    dx_list(k) = dx;
    

    % Define information about the step sizes in the "time" dimension
    dt = dx^2 * ratio(k);
    T = 0.01;
    N = floor(T * (1.0 + 1.0e-10) / dt);

    alpha = 1;

    % The initial condition
    psi_n = 4*sin(4*pi*X).*sin(5*pi*Y)+sin(pi*X).*sin(pi*Y);

    psi_np1 = zeros(size(psi_n));

    er_norm1 = 0;
    er_norm2 = 0;

    for n = 0:N - 1
    psi_np1(2:end-1,2:end-1) = psi_n(2:end-1,2:end-1) - dt*psi_n(2:end-1,2:end-1)...
        + (dt/(dx^2))*(1-alpha)*(psi_n(2:end-1,3:end)+psi_n(1:end-2,2:end-1)...
        - 4*psi_n(2:end-1,2:end-1) + psi_n(3:end,2:end-1)+ psi_n(2:end-1,1:end-2))...
        + (dt/(2*dx^2))*alpha*(psi_n(1:end-2,3:end)+psi_n(3:end,3:end)-4*psi_n(2:end-1,2:end-1)...
        + psi_n(1:end-2,1:end-2) + psi_n(3:end,1:end-2));

  
    % Exact solution after (n+1) steps
    Psi_np1 = 4*exp(-dt*(n+1)-(16+25)*pi^2*(n+1)*dt).*sin(4*pi*X).*sin(5*pi*Y) ...
        + exp(-dt*(n+1)-(1+1)*(n+1)*pi^2*dt).*sin(pi*X).*sin(pi*Y);
  
    err = (psi_np1 - Psi_np1);
    
    er_norm1 = max(er_norm1 ,sqrt(sum(sum(err.^2))));
    er_norm2 = max(er_norm2, sqrt(dx*dx.*sum(sum(err.^2))));
    
    psi_n = psi_np1;
    end
    err_norm1_list(k) = abs(max(er_norm1));
    err_norm2_list(k) = abs(max(er_norm2));

end
figure(1)
plot(ratio, err_norm2_list, '-o')
xlabel({'$\Delta t/\Delta x^2$'}, 'Interpreter', 'latex')
ylabel('Error')
title('Error Vs \Delta t/\Delta x^2')  