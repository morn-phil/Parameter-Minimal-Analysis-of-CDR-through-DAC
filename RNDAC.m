function myRNDAC()
    % Clear workspace
    clear; clc; close all;
    
    % --- Define parameters ---
    % Rate constants
    k1 = 0.4; k2 = 0.3; k4 = 0.1; k5 = 0.05; k6 = 0.01;
    am = 0.05; b = 1;
    
    % Reaction orders
    p1 = -0.2; q1 = 4;  
    p2 = 3; q2 = -0.5;  
    
    % Pack all parameters into a structure
    params.k1 = k1; params.k2 = k2; params.k4 = k4; params.k5 = k5; params.k6 = k6;
    params.am = am; params.b = b;
    params.p1 = p1; params.q1 = q1; params.p2 = p2; params.q2 = q2;
    
    % --- Time span for simulation ---
    tspan = [0 100];
    
    % --- Multiple initial conditions ---
    initial_conditions = [
        0.25, 0.25, 0.25, 0.25, 0;  % IC 1
        0.35, 0.15, 0.35, 0.15, 0;  % IC 2  
        0.15, 0.05, 0.45, 0.35, 0   % IC 3
    ];
    
    % --- Colors for different ICs ---
    colors = [0.8500, 0.3250, 0.0980;  % Orange
              0.0000, 0.4470, 0.7410;  % Blue
              0.4660, 0.6740, 0.1880]; % Green
    
    % --- Create figure with subplots ---
    figure('Position', [100, 100, 1400, 800]);
    
    % Variable names
    var_names = {'A_1', 'A_2', 'A_3', 'A_4', 'A_5'};
    
    % Loop through initial conditions
    for ic = 1:size(initial_conditions, 1)
        y0 = initial_conditions(ic, :)';
        
        % Solve ODE system
        [t, y] = ode45(@(t, y) chemicalODEs(t, y, params), tspan, y0);
        
        % Plot each variable
        for var = 1:5
            subplot(2, 3, var);
            hold on;
            plot(t, y(:, var), 'Color', colors(ic, :), 'LineWidth', 2, ...
                 'DisplayName', ['IC ', num2str(ic)]);
            title(['Variable: ', var_names{var}]);
            xlabel('Time');
            ylabel(['Concentration ', var_names{var}]);
            grid on;
            legend('Location', 'best');
        end
    end
    
    % Add phase plot
    % subplot(2, 3, 6);
    % hold on;
    % for ic = 1:size(initial_conditions, 1)
    %     y0 = initial_conditions(ic, :)';
    %     [t, y] = ode45(@(t, y) chemicalODEs(t, y, params), tspan, y0);
    %     plot(y(:, 1), y(:, 2), 'Color', colors(ic, :), 'LineWidth', 2, ...
    %          'DisplayName', ['IC ', num2str(ic)]);
    % end
    % title('Phase Plot: x_1 vs x_2');
    % xlabel('x_1');
    % ylabel('x_2');
    % grid on;
    % legend('Location', 'best');
    
    sgtitle('RNDAC');
end

% --- ODE Function Definition ---
function dxdt = chemicalODEs(t, x, params)
    % Extract parameters
    k1 = params.k1; k2 = params.k2; k4 = params.k4; 
    k5 = params.k5; k6 = params.k6; am = params.am; b = params.b;
    p1 = params.p1; q1 = params.q1; p2 = params.p2; q2 = params.q2;
    
    % Initialize derivatives
    dxdt = zeros(5, 1);
    
    % ODE equations 
    dxdt(1) = k1 * x(1)^(p1) * x(2)^(q1) - k2 * x(1)^(p2) * x(2)^(q2);
    dxdt(2) = - k1*x(1)^(p1)*x(2)^(q1) + k2*x(1)^(p2)*x(2)^(q2) - am*x(2) + am*b*x(3) + k4*x(4) -k5*x(2);
    dxdt(3) = am * x(2) - am * b * x(3);
    dxdt(4) = k6 * x(5) - k4 * x(4);
    dxdt(5) = k5 * x(2) - k6 * x(5);
end