clear all;
close all;
clc;
timerVal = tic;
startTime = cputime;

%% Define parameters
H = 0.4; % Height of the plate (m)
L = 0.3; % Length of the plate (m)
nx = 31; % Number of grid points in x-direction
ny = 41; % Number of grid points in y-direction
dx = L / (nx - 1); % Grid spacing in x-direction
dy = H / (ny - 1); % Grid spacing in y-direction
alpha = 11.234e-5; % Thermal diffusivity (m^2/s)
dt = 0.1; % Time step size (s)
gamma_x = alpha * dt / dx^2; % Stability parameter in x-direction
gamma_y = alpha * dt / dy^2; % Stability parameter in y-direction

% Boundary conditions
T1 = 40; % Bottom boundary temperature (°C)
T2 = 0;  % Right boundary temperature (°C)
T3 = 10; % Top boundary temperature (°C)
T4 = 0;  % Left boundary temperature (°C)

%% Initialize the temperature field
T = zeros(ny, nx); % Initial temperature field

% Apply boundary conditions
T(1, :) = T1;       % Bottom boundary
T(end, :) = T3;     % Top boundary
T(:, end) = T2;     % Right boundary
T(:, 1) = T4;       % Left boundary

%% Convergence criteria and error tracking
error_mag = 1;
error_req = 0.01;
iterations = 0;

% Store temperature values at specific times
time_steps_to_store = [0, 0.1, 0.5, 1, 5, 10, 20, 40];
time_steps_to_store_indices = round(time_steps_to_store / dt) + 1;
stored_temperatures_x15_y_all = zeros(length(time_steps_to_store), ny);
stored_temperatures_y20_x_all = zeros(length(time_steps_to_store), nx);

%% ADI Method Implementation with TDMA Solver
while true
    T_old = T;

    % Step 1: Solve in x-direction (Implicit in x, Explicit in y)
    for j = 2:ny-1
        a_x = -gamma_x * ones(1, nx-2);        % Sub-diagonal coefficients
        b_x = (1 + 2*gamma_x) * ones(1, nx-2); % Main diagonal coefficients
        c_x = -gamma_x * ones(1, nx-2);        % Super-diagonal coefficients

        d_x = T(j, 2:end-1)' + gamma_y * (T(j+1, 2:end-1)' + T(j-1, 2:end-1)' - 2*T(j, 2:end-1)');

        d_x(1) = d_x(1) + gamma_x * T(j, 1);       % Left boundary contribution
        d_x(end) = d_x(end) + gamma_x * T(j, end); % Right boundary contribution

        T(j, 2:end-1) = TDMA(a_x, b_x, c_x, d_x);  % Solve using TDMA
    end

    % Step 2: Solve in y-direction (Implicit in y, Explicit in x)
    for i = 2:nx-1
        a_y = -gamma_y * ones(1, ny-2);        % Sub-diagonal coefficients
        b_y = (1 + 2*gamma_y) * ones(1, ny-2); % Main diagonal coefficients
        c_y = -gamma_y * ones(1, ny-2);        % Super-diagonal coefficients

        d_y = T(2:end-1, i) + gamma_x * (T(2:end-1, i+1) + T(2:end-1, i-1) - 2*T(2:end-1, i));

        d_y(1) = d_y(1) + gamma_y * T(1, i);       % Bottom boundary contribution
        d_y(end) = d_y(end) + gamma_y * T(end, i); % Top boundary contribution

        T(2:end-1, i) = TDMA(a_y', b_y', c_y', d_y'); % Solve using TDMA
    end

    iterations = iterations + 1;

    %% Convergence Check and Data Storage
    error_mag = max(max(abs(T - T_old)));

    if ismember(iterations, time_steps_to_store_indices)
        index_to_store = find(time_steps_to_store_indices == iterations);
        stored_temperatures_x15_y_all(index_to_store, :) = T(:, round(nx/2)); % At x=0.15 m along y-direction
        stored_temperatures_y20_x_all(index_to_store, :) = T(round(ny/2), :); % At y=0.2 m along x-direction
    end

    if error_mag < error_req && iterations > max(time_steps_to_store_indices)
        break;
    end
end

%% Plot Steady-State Temperature Distribution
x_dom = linspace(0, L, nx);
y_dom = linspace(0, H, ny);
[X, Y] = meshgrid(x_dom, flip(y_dom));

figure;
contourf(X, Y, flipud(T), 'LineColor', 'none');
colorbar;
xlabel('x (m)');
ylabel('y (m)');
title('Steady-State Temperature Distribution');

%% Plot Temperature Profiles at Specific Times
% Plot temperature along y-direction at x=0.15 m for different times
figure;
for k = 1:length(time_steps_to_store)
    plot(flip(y_dom), stored_temperatures_x15_y_all(k,:), 'DisplayName', ['t=' num2str(time_steps_to_store(k)) 's']);
    hold on;
end
xlabel('y (m)');
ylabel('Temperature (°C)');
title('Temperature Profiles Along y-Direction at x=0.15 m');
legend show;

% Plot temperature along x-direction at y=0.2 m for different times
figure;
for k = 1:length(time_steps_to_store)
    plot(x_dom, stored_temperatures_y20_x_all(k,:), 'DisplayName', ['t=' num2str(time_steps_to_store(k)) 's']);
    hold on;
end
xlabel('x (m)');
ylabel('Temperature (°C)');
title('Temperature Profiles Along x-Direction at y=0.2 m');
legend show;

Total_computational_Time = toc(timerVal);
disp(Total_computational_Time)
CPU_time_1 = cputime - startTime;
disp(CPU_time_1)

%% TDMA Solver Function (local functions must be at the end of a script file)
function x = TDMA(a, b, c, d)
    n = length(d);
    c_star = zeros(1, n);
    d_star = zeros(1, n);
    x = zeros(1, n);

    c_star(1) = c(1) / b(1);
    d_star(1) = d(1) / b(1);

    for i = 2:n
        temp = b(i) - a(i) * c_star(i-1);
        c_star(i) = c(i) / temp;
        d_star(i) = (d(i) - a(i) * d_star(i-1)) / temp;
    end

    x(n) = d_star(n);
    for i = n-1:-1:1
        x(i) = d_star(i) - c_star(i) * x(i+1);
    end
end
