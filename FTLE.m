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
dt = 0.001; % Time step size (s)
gamma_x = alpha * dt / dx^2; % Stability parameter in x-direction
gamma_y = alpha * dt / dy^2; % Stability parameter in y-direction

% Boundary conditions
T1 = 40; % Bottom boundary temperature (°C)
T2 = 0;  % Right boundary temperature (°C)
T3 = 10; % Top boundary temperature (°C)
T4 = 0;  % Left boundary temperature (°C)

%% Initialize the temperature field
T = zeros(ny, nx); % Initial temperature field
T_new = T;

% Apply boundary conditions
T(1, :) = T1;       % Bottom boundary
T(end, :) = T3;     % Top boundary
T(:, end) = T2;     % Right boundary
T(:, 1) = T4;       % Left boundary

T_new = T;

%% Convergence criteria and error tracking
error_mag = 1;
error_req = 0.01;
iterations = 0;

% Store temperature values at specific times
time_steps_to_store = [0, 0.1, 0.5, 1, 5, 10, 20, 40];
time_steps_to_store_indices = round(time_steps_to_store / dt) + 1;
stored_temperatures_x15_y_all = zeros(length(time_steps_to_store), ny);
stored_temperatures_y20_x_all = zeros(length(time_steps_to_store), nx);

%% Time-stepping loop for unsteady diffusion equation
while true
    for i = 2:ny-1
        for j = 2:nx-1
            T_new(i, j) = T(i, j) + gamma_y * (T(i, j+1) + T(i, j-1) - 2*T(i, j)) ...
                                    + gamma_x * (T(i+1, j) + T(i-1, j) - 2*T(i, j));
        end
    end

    iterations = iterations + 1;

    % Calculate error magnitude for convergence check
    error_mag = max(max(abs(T_new - T)));
    if error_mag < error_req && iterations > max(time_steps_to_store_indices)
        break;
    end

    % Store temperatures at specific time steps for plotting later
    if ismember(iterations, time_steps_to_store_indices)
        index_to_store = find(time_steps_to_store_indices == iterations);
        stored_temperatures_x15_y_all(index_to_store, :) = T_new(:, round(nx/2)); % At x=0.15 m along y-direction
        stored_temperatures_y20_x_all(index_to_store, :) = T_new(round(ny/2), :); % At y=0.2 m along x-direction
    end

    % Update old temperature field with new values
    T = T_new;
end

%% Plot steady-state temperature distribution
x_dom = linspace(0, L, nx);
y_dom = linspace(0, H, ny);
[X, Y] = meshgrid(x_dom, flip(y_dom)); % Flip y-axis for proper orientation

figure;
contourf(X, Y, flipud(T), 12);
colorbar;
xlabel('x (m)');
ylabel('y (m)');
title('Steady-State Temperature Distribution at time t');

%% Plot temperature profiles at specific times
% Plot temperature along y-direction at x=0.15 m for different times
figure;
for k = 1:length(time_steps_to_store)
    plot(flip(y_dom), stored_temperatures_x15_y_all(k, :), 'DisplayName', ['t=' num2str(time_steps_to_store(k)) 's']);
    hold on;
end
xlabel('y (m)');
ylabel('Temperature (°C)');
title('Temperature Profiles Along y-Direction at x=0.15 m');
legend show;

% Plot temperature along x-direction at y=0.2 m for different times
figure;
for k = 1:length(time_steps_to_store)
    plot(x_dom, stored_temperatures_y20_x_all(k, :), 'DisplayName', ['t=' num2str(time_steps_to_store(k)) 's']);
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
