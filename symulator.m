function symulator
    
end

function [t, y] = metoda_ulepszona_eulera(f, t0, tn, y0, h)
    t = t0:h:tn;
    N = length(t);

    y = zeros(length(y0), N);
    y(:,1) = y0;

    for n = 1:(N-1)
        y(:, n + 1) = y(:, n) + h * f(t(n) + h/2, y(:,n) + h/2 * f(t(n), y(:, n)));
    end
end

function [t, y] = metoda_prosta_eulera(f, t0, tn, y0, h)
    t = t0:h:tn;
    N = length(t);
    
    y = zeros(length(y0), N);
    y(:,1) = y0;

    for n = 1:(N-1)
        y(:, n + 1) = y(:, n) + h * f(t(n), y(:,n));
    end
end