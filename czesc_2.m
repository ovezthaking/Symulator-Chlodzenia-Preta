function czesc_2
    format long;
    clf;
    dT_tab = [-1500, -1000, -300, -50, -1, 1,20, 50, 200, 400, 1000, 2000];
    h_tab  = [178, 176, 168, 161, 160, 160, 160.2, 161, 165, 168, 174, 179];
    
    hold on;
    plot(dT_tab, h_tab, 'o')
    dT_plot = linspace(dT_tab(1), dT_tab(end), 2000);
    for m=2:8
        h_approx = approx(dT_tab, h_tab, m, dT_plot);

        plot(dT_plot, h_approx, "DisplayName", sprintf("m = %d", m))
    end
    ylim([min(h_tab)-1, max(h_tab)+1])
    legend("Location","southwest")
    hold off;
    

    figure;
    h_interp = zeros(length(dT_plot), 1);

    for i = 1:length(dT_plot)
        h_interp(i) = interp_newton(dT_tab, h_tab, dT_plot(i));
    end
    
    [dT_eq, h_eq] = generuj_wezly_rownoodlegle(dT_tab, h_tab, 50);
    h_spline = spline(dT_eq, h_eq, dT_plot);

    h_approx = approx(dT_tab, h_tab, 6, dT_plot);

    plot(dT_tab, h_tab, 'o')
    hold on;
    plot(dT_plot, h_interp, 'k', 'LineWidth', 1.5)
    plot(dT_plot, h_spline)
    plot(dT_plot, h_approx)
    grid on;
    legend('Dane', 'Interpolacja Newtona', 'Location', 'southwest')
    ylim([min(h_tab)-1, max(h_tab)+1])
    hold off;
    
end

function [x_eq, y_eq] = generuj_wezly_rownoodlegle(x_tab, y_tab, hx)
    x_eq = x_tab(1):hx:x_tab(end);
    y_eq = zeros(size(x_eq));

    for k = 1:length(x_eq)
        xk = x_eq(k);

        i = 1;
        while i < length(x_tab) - 1 && xk > x_tab(i + 1)
            i = i + 1;
        end

        y_eq(k) = y_tab(i) + (y_tab(i + 1) - y_tab(i)) * ...
                  (xk - x_tab(i)) / (x_tab(i + 1) - x_tab(i));
    end
end