%% initial setting
[mu, ~, ~, ~]      = parameter(2); % Earth-Moon
[L1, L2, L3, L4, L5] = librationPoints(mu);

count_max     = 8000;
plot_interval = 100;

Jacobi_min = 3.0;
Jacobi_max = 3.1723;

Jacobi       = zeros(count_max/plot_interval, 1);
x0_corrected = zeros(count_max/plot_interval, 6);
t0_corrected = zeros(count_max/plot_interval, 1);

C_xn = 3.12;
% C_xn = Jacobi_const([L2;zeros(3,1)],mu);
C_error = 1e-13;
disp(C_xn);

%continuation-----------------------------------------------------------------------------------
iteration_max = 100;
threshold     = 1e-10;
delta         = 4e-5;
count         = 0;
count_max     = 8000;
detect_period_orbit_1 = zeros(8,1);

x0_1 = [0.836900082907655, 0, 0, 0, 1.770874936727959e-06, 0];
t0_1 = 1.3458;

%Loop for first orbit------------------------------------------------------------------------------------------------
while 1
    count = count + 1;

    % Differential correction-----------------------------------------------
    for iteration = 1:iteration_max
        [x_n_1, t_n_1, C] = fun_differential_correction_cr3bp(x0_1, t0_1, mu, options_ODE);

        tspan = [0 2*t_n_1];
        [t_corrected_1, x_corrected_1] = ode113(@(t, x) fun_cr3bp(t, x, mu), tspan, x_n_1, options_ODE);

        x_error_1 = norm(x_corrected_1(end, :) - x_corrected_1(1, :));

        if x_error_1 < threshold
            break;
        end

        if x_error_1 > 1e+3
            disp('calculation diverged');
            return;
        end

        if iteration == iteration_max
            count = count - 1;
            x0_1(1) = x0_1(1) - delta;
            delta = delta/10;
        end

        x0_1 = x_n_1;
        t0_1 = t_n_1;
    end

    if C < C_xn
        x0_1(1) = x0_1(1) - delta;
        delta = delta/5;
        disp(strcat('delta changed : count = ', num2str(count)));

    end

    if(count >= 100)&&(mod(count,100)==0)
        disp(strcat('count = ', num2str(count)));
    end

    if abs(C-C_xn)<C_error
        additional = [C; t_n_1; x_n_1(:,1)];
        detect_period_orbit_1 = additional;
        break;
    end

    x0_1(1) = x0_1(1) + delta;

    if count == count_max - 1
        break
    end
end

disp('L1-sucucessfuly finished');