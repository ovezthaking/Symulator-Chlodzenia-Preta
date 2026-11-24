function czesc_1
    format long;
    clc; close all;

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
    legend('Location', 'northeast', 'FontSize', 12);
    xlim([0, 5]); ylim([0, 1250]);
    hold off;

    przypadki = [
        1200, 25,  2.5;
        800, 25,  2.5;
        1100, 70,  2.5;
        1100, 70,  5.0;
        1100, 70, 10.0
    ];
    
    czasy_odczytu = [2, 3, 4, 5];
    dt_sim = 0.001;

    fprintf('\n--- WYNIKI SYMULACJI DLA WARIANTÓW Z EKSPERYMENTU ---\n');
    for i = 1:size(przypadki, 1)
        Tb0_i = przypadki(i, 1);
        Tw0_i = przypadki(i, 2);
        mw_i  = przypadki(i, 3);

        [ti, Tbi, Twi] = symulator(h, A, mb, mw_i, cb, cw, ...
            [Tb0_i; Tw0_i], [0, 5], dt_sim);

        fprintf('\nWariant %d (Tb(0) = %d st.C, Tw(0) = %d st.C, mw = %.1f kg):\n', ...
            i, Tb0_i, Tw0_i, mw_i);

        for t_odcz = czasy_odczytu
            idx = round(t_odcz / dt_sim) + 1;
            fprintf('  t = %d s  ->  Tb = %7.2f st.C,  Tw = %7.2f st.C\n', ...
                t_odcz, Tbi(idx), Twi(idx));
        end
    end


    figure(2);
    kolory = lines(size(przypadki, 1));

    subplot(1, 2, 1);
    hold on; box on; grid on;
    for i = 1:size(przypadki, 1)
        Tb0_i = przypadki(i, 1);
        Tw0_i = przypadki(i, 2);
        mw_i  = przypadki(i, 3);

        [ti, Tbi, ~] = symulator(h, A, mb, mw_i, cb, cw, [Tb0_i; Tw0_i], [0, 5], 0.01);
        plot(ti, Tbi, 'Color', kolory(i, :), 'LineWidth', 1.2, 'DisplayName', ...
            sprintf('T_{b0}=%d, T_{w0}=%d, m_w=%.1f kg', Tb0_i, Tw0_i, mw_i));
    end
    xlabel('t [s]'); ylabel('T_b [^\circC]');
    title('Temperatura pręta T_b(t)');
    legend('Location', 'northeast', 'FontSize', 10);
    hold off;

    subplot(1, 2, 2);
    hold on; box on; grid on;
    for i = 1:size(przypadki, 1)
        Tb0_i = przypadki(i, 1);
        Tw0_i = przypadki(i, 2);
        mw_i  = przypadki(i, 3);

        [ti, ~, Twi] = symulator(h, A, mb, mw_i, cb, cw, [Tb0_i; Tw0_i], [0, 5], 0.01);
        plot(ti, Twi, 'Color', kolory(i, :), 'LineWidth', 1.2, 'DisplayName', ...
            sprintf('T_{b0}=%d, T_{w0}=%d, m_w=%.1f kg', Tb0_i, Tw0_i, mw_i));
    end
    xlabel('t [s]'); ylabel('T_w [^\circC]');
    title('Temperatura oleju T_w(t)');
    legend('Location', 'southeast', 'FontSize', 10);
    hold off;


    figure('Name', 'Czesc 1 - Wplyw kroku dt', 'Position', [150, 150, 950, 420]);
    kroki_test = [0.2, 0.1, 0.05, 0.001];

    subplot(1, 2, 1);
    hold on; box on; grid on;
    for dt_test = kroki_test
        [t_dt, Tb_dt, ~] = symulator(h, A, mb, mw, cb, cw, T0, [0, 2.5], dt_test, 'euler');
        plot(t_dt, Tb_dt, 'LineWidth', 1.2, 'DisplayName', sprintf('dt = %.3f s', dt_test));
    end
    plot(t_45, T_45(:, 1), 'k--', 'LineWidth', 1.5, 'DisplayName', 'ode45');
    xlim([0, 2.5]); ylim([0, 1200]);
    xlabel('t [s]'); ylabel('T_b [^\circC]');
    title('Metoda prosta Eulera O(\Deltat)');
    legend('Location', 'northeast', 'FontSize', 9);
    hold off;

    subplot(1, 2, 2);
    hold on; box on; grid on;
    for dt_test = kroki_test
        [t_dt, Tb_dt, ~] = symulator(h, A, mb, mw, cb, cw, T0, [0, 2.5], dt_test, 'euler_improved');
        plot(t_dt, Tb_dt, 'LineWidth', 1.2, 'DisplayName', sprintf('dt = %.3f s', dt_test));
    end
    plot(t_45, T_45(:, 1), 'k--', 'LineWidth', 1.5, 'DisplayName', 'ode45');
    xlim([0, 2.5]); ylim([0, 1200]);
    xlabel('t [s]'); ylabel('T_b [^\circC]');
    title('Metoda ulepszona Eulera O(\Deltat^2)');
    legend('Location', 'northeast', 'FontSize', 9);
    hold off;


    fprintf('\n--- ANALIZA DOKŁADNOŚCI METOD W ZALEŻNOŚCI OD KROKU dt ---\n');
    kroki_err = [0.1, 0.05, 0.02, 0.01, 0.001];
    for dt_e = kroki_err
        [te, Tb_e, ~]  = symulator(h, A, mb, mw, cb, cw, T0, interval, dt_e, ...
            'euler');
        [tu, Tb_u, ~]  = symulator(h, A, mb, mw, cb, cw, T0, interval, dt_e, ...
            'euler_improved');

        Tb_ref = interp1(t_45, T_45(:, 1), te, 'spline');

        err_e = max(abs(Tb_e - Tb_ref));
        err_u = max(abs(Tb_u - Tb_ref));
        fprintf(['dt = %.3f s -> Max błąd Euler prosty = %7.4f st.C | ' ...
            'Euler ulepszony = %7.4f st.C\n'], ...
            dt_e, err_e, err_u);
    end
end