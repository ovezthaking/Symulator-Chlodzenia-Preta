function czesc_1
    format long;

    h = 160;
    A = 0.0109;
    mb = 0.2;
    mw = 2.5;
    cb = 3.85;
    cw = 4.1813;

    T0 = [1200; 25];
    interval = [0, 5];
    dt = 0.05;

    [t_mpe, Tb_mpe, Tw_mpe] = symulator(h, A, mb, mw, cb, cw, T0, ...
        interval, dt, "euler");
    [t_mpu, Tb_mpu, Tw_mpu] = symulator(h, A, mb, mw, cb, cw, T0, ...
        interval, dt);

    F_ode = @(~, T) [-(h * A) / (mb * cb) * (T(1) - T(2)); ...
        (h * A) / (mw * cw) * (T(1) - T(2))];

    [t_23, T_23] = ode23(F_ode, interval, T0);
    opts_45 = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    [t_45, T_45] = ode45(F_ode, interval, T0, opts_45);

    figure('Name', 'Czesc 1 - Przebieg bazowy', 'Position', [100, 100, 850, 500]);
    hold on; box on; grid on;
    plot(t_mpu, Tb_mpu, '-x', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.2, ...
        'MarkerSize', 6, 'DisplayName', 'T_b(t) - Zmodyfikowany Euler (dt = 0.05 s)');
    plot(t_mpu, Tw_mpu, '-o', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.2, ...
        'MarkerSize', 6, 'DisplayName', 'T_w(t) - Zmodyfikowany Euler (dt = 0.05 s)');
    plot(t_mpe, Tb_mpe, 'b--', 'LineWidth', 1.1, 'DisplayName', 'T_b(t) - Prosty Euler (dt = 0.05 s)');
    plot(t_23, T_23(:, 1), 'g:', 'LineWidth', 1.4, 'DisplayName', 'T_b(t) - Wbudowany ode23');
    plot(t_45, T_45(:, 1), 'k:', 'LineWidth', 1.5, 'DisplayName', 'T_b(t) - Wbudowany ode45');
    xlabel('t [s]', 'FontSize', 11);
    ylabel('T [^\circC]', 'FontSize', 11);
    title('Przebieg czasowy temperatury pręta T_b(t) i oleju T_w(t)', 'FontSize', 12);
    legend('Location', 'northeast', 'FontSize', 9);
    xlim([0, 5]); ylim([0, 1250]);
    hold off;
end