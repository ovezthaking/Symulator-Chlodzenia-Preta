function [t, Tb, Tw] = symulator(h, A, mb, mw, cb, cw, T0, interval, dt, metoda)
    % symulator - Główny symulator procesu chłodzenia pręta w oleju
    % Wejście:
    %   h         - współczynnik przewodnictwa cieplnego
    %   A         - powierzchnia pręta
    %   mb        - masa pręta
    %   mw        - masa oleju chłodzącego
    %   cb        - pojemność cieplna metalu
    %   cw(T)     - pojemność cieplna oleju chłodzącego zależna
    %               od temperatury T (nieliniowe)
    %   T0        - warunek początkowy [Tb(0); Tw(0)]
    %   interval  - przedział czasu [t_start, t_end]
    %   dt        - krok całkowania [s]
    %   metoda    - 'euler' lub 'euler_modified' 
    %               (domyślnie 'euler_modified')
    % Wyjście:
    %   t         - wektor chwil czasowych [s]
    %   Tb        - przebieg temperatury pręta [st. C]
    %   Tw        - przebieg temperatury oleju [st. C]

    if nargin < 10
        metoda = "euler_improved";
    end
    
    F = @(T) rownania_stanu(T, h, A, mb, mw, cb, cw);

    switch lower(metoda)
        case 'euler'
            [t, Y] = metoda_prosta_eulera(F, interval(1), ...
                interval(end), T0, dt);
        case 'euler_improved'
            [t, Y] = metoda_ulepszona_eulera(F, interval(1), ...
                interval(end), T0, dt)
        otherwise
            error("Nieznana metoda - wybierz 'euler' lub" + ...
                "'euler_improved'")
    end

    Tb = Y(1, :);
    Tw = Y(2, :);
end

function dT = rownania_stanu(T, h, A, mb, mw, cb, cw)
     Tb = T(1);
     Tw = T(2);

     dTb_dt = - (h*A/mb*cb) * (Tb-Tw);
     dTw_dw = h*A/mw*cw * (Tb-Tw);

     dT = [dTb_dt, dTw_dw];
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