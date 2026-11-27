function czesc_2
    format long;
    clf;
    dT_tab = [-1500, -1000, -300, -50, -1, 1,20, 50, 200, 400, 1000, 2000];
    h_tab  = [178, 176, 168, 161, 160, 160, 160.2, 161, 165, 168, 174, 179];
    
    % Aproksymacja - dobór stopnia wielomianu

    fprintf('DOBOR STOPNIA WIELOMIANU APROKSYMACJI\n');
    for m = 1:8
        h_approx = approx(dT_tab, h_tab, m, dT_tab);
        err_abs = abs(h_approx - h_tab);
        bl = sqrt(mean((h_tab - h_approx).^2));
        max_rel = max(err_abs ./ h_tab) * 100;
        fprintf('Stopien %d: bl = %7.4f | MaxErr = %7.4f (%.3f%%)\n', ...
            m, bl, max(err_abs), max_rel);
    end
    
    % Wykres 1
    figure(1);
    hold on;
    grid on;
    box on;

    plot(dT_tab, h_tab, 'o')
    dT_plot = linspace(dT_tab(1), dT_tab(end), 2000);
    for m=2:8
        h_approx = approx(dT_tab, h_tab, m, dT_plot);

        plot(dT_plot, h_approx, "DisplayName", sprintf("m = %d", m));
    end
    
    % xlim([0, 1200]);
    ylim([min(h_tab)-1, 185]);
    xlabel('\DeltaT [^\circC]'); ylabel('h [W\cdotm^{-2}]');
    title(['Aproksymacja metodą najmniejszych kwadratów' ...
        ' dla roznych stopni wielomianu']);
    legend('Location', 'southwest', 'FontSize' ,15);
    hold off;
    

    % Dobór metody
    m_optymalne = 5;

    [dT_rown, h_rown] = generuj_wezly_rownoodlegle(dT_tab, h_tab, 50);
    
    h_interp = interp_newton(dT_tab, h_tab, dT_tab);
    h_approx = approx(dT_tab, h_tab, m_optymalne, dT_tab);
    h_spline = splajn(dT_rown, h_rown, dT_tab);
    
    fprintf("\nPOROWNANIE METOD NUMERYCZNYCH\n")
    metody = [h_interp; h_approx; h_spline];
    nazwy = ["Interpolacja", "Aproksymacja", "Splajn"];

    for i=1:size(metody, 1)
        err_abs = abs(h_tab - metody(i,:));
        bl = sqrt(mean((h_tab - metody(i,:)).^2));
        max_rel = max(err_abs ./ h_tab) * 100;
        fprintf("%s: bl = %.4e | MaxErr = %.4e (%.4f%%)\n", ...
            nazwy(i), bl, max(err_abs), max_rel);
    end


    % Wykres 2
    % Interpolacja
    h_interp_plot = interp_newton(dT_tab, h_tab, dT_plot);
    
    % Splajn
    h_spline_plot = splajn(dT_rown, h_rown, dT_plot);
    
    % Aproksymacja
    h_approx_plot = approx(dT_tab, h_tab, m_optymalne, dT_plot);
    
    figure(2);
    plot(dT_tab, h_tab, 'o', "DisplayName", "Dane");
    hold on;
    grid on;
    box on;

    plot(dT_plot, h_interp_plot, 'k', 'LineWidth', 1.5, ...
        "DisplayName", "Interpolacja Newtona");
    plot(dT_plot, h_spline_plot, 'LineWidth', 1, ...
        "DisplayName", "Funkcja sklejana");
    plot(dT_plot, h_approx_plot, 'LineWidth', 1.5, ...
        "DisplayName", "Aproksymacja mnk");
    
    xlabel('\DeltaT [^\circC]'); ylabel('h [W\cdotm^{-2}]');
    title(['Metody numeryczne w pelnym zakresie' ...
        ' \DeltaT \in [-1500, 2000]^\circC']);
    legend('Location', 'southwest', 'FontSize' ,15);
    % xlim([0, 1200]);
    ylim([min(h_tab)-1, 185]);
    hold off;

    % Wykres 3
    h_approx_plot_m4 = approx(dT_tab, h_tab, 4, dT_plot);
    
    figure(3);
    plot(dT_tab, h_tab, 'o', "DisplayName", "Dane pomiarowe");
    hold on;
    grid on;
    box on;

    plot(dT_plot, h_interp_plot, 'k', 'LineWidth', 1.5, ...
        "DisplayName", "Interpolacja Newtona");
    plot(dT_plot, h_spline_plot, 'LineWidth', 1, ...
        "DisplayName", "Funkcja sklejana");
    plot(dT_plot, h_approx_plot, 'LineWidth', 1.5, ...
        "DisplayName", "Aproksymacja mnk m=5");
    plot(dT_plot, h_approx_plot_m4, 'LineWidth', 1.5, ...
        "DisplayName", "Aproksymacja mnk m=4");

    xlabel('\DeltaT [^\circC]'); ylabel('h [W\cdotm^{-2}]');
    title(['Metody numeryczne zakresie roboczym symulacji' ...
        ' \DeltaT \in [0, 1200]^\circC']);
    legend('Location', 'northwest', 'FontSize' ,15);
    xlim([0, 1200]);
    ylim([min(h_tab)-1, 185]);
    hold off;
end

function [x, y] = generuj_wezly_rownoodlegle(x_tab, y_tab, hx)
    x = x_tab(1):hx:x_tab(end);
    y = zeros(size(x));

    for k = 1:length(x)
        xk = x(k);

        i = 1;
        while i < length(x_tab) - 1 && xk > x_tab(i + 1)
            i = i + 1;
        end

        y(k) = y_tab(i) + (y_tab(i + 1) - y_tab(i)) * ...
                  (xk - x_tab(i)) / (x_tab(i + 1) - x_tab(i));
    end
end
