function Sx = spline(x, y, x_target)
    [C, a, b, h] = wyznacz_wspolczynniki_splajnu(x, y);

    Sx = zeros(size(x_target));
    for k = 1:numel(x_target)
        xk = min(max(x_target(k), a), b);
        phi = oblicz_tablice_phi(C, a, h, xk);
        Sx(k) = oblicz_wartosc_splajnu(phi, C);
    end
end

function phi = oblicz_tablice_phi(C, a, h, xValue)
    phi = zeros(length(C), 1);
    x = @(i) round(a + i * h, 10);
    invh3 = 1 / (h^3);

    for i = -1:length(C) - 2
        if x(i - 2) <= xValue && xValue <= x(i - 1)
            phi(i + 2) = invh3 * ((xValue - x(i - 2))^3);
        elseif x(i - 1) <= xValue && xValue <= x(i)
            phi(i + 2) = invh3 * ((xValue - x(i - 2))^3 - 4 * (xValue - x(i - 1))^3);
        elseif x(i) <= xValue && xValue <= x(i + 1)
            phi(i + 2) = invh3 * ((x(i + 2) - xValue)^3 - 4 * (x(i + 1) - xValue)^3);
        elseif x(i + 1) <= xValue && xValue <= x(i + 2)
            phi(i + 2) = invh3 * ((x(i + 2) - xValue)^3);
        end
    end
end

function Sx = oblicz_wartosc_splajnu(phi, C)
    Sx = 0;
    for i = 1:length(phi)
        Sx = Sx + phi(i) * C(i);
    end
end

function [C, a, b, h] = wyznacz_wspolczynniki_splajnu(xi, yi)
    n = length(xi) - 1;
    a = xi(1);
    b = xi(end);
    h = (b - a) / n;

    alpha = (yi(2) - yi(1)) / h;
    beta  = (yi(end) - yi(end-1)) / h;

    d = oblicz_tablice_d(xi, yi, alpha, beta, h);
    A = znajdz_macierz_A(d);
    C = oblicz_tablice_wspolczynnikow(A, d, h, n, alpha, beta);
end

function d = oblicz_tablice_d(x, y, alpha, beta, h)
    if length(x) ~= length(y)
        error("Tablica x oraz y muszą być tej samej wielkości");
    end
    d = y(:);
    d(1)   = y(1) + (h / 3) * alpha;
    d(end) = y(end) - (h / 3) * beta;
end

function A = znajdz_macierz_A(d)
    A = zeros(length(d), length(d));
    A(1, 1) = 4;
    A(1, 2) = 2;
    A(end, end) = 4;
    A(end, end - 1) = 2;
    for i = 2:length(A) - 1
        A(i, i - 1) = 1;
        A(i, i)     = 4;
        A(i, i + 1) = 1;
    end
end

function C = oblicz_tablice_wspolczynnikow(A, d, h, n, alpha, beta)
    C = A \ d(:);
    c = @(i) C(i + 1);
    c_m1  = c(1) - (h / 3) * alpha;
    c_np1 = c(n - 1) + (h / 3) * beta;
    C = [c_m1; C; c_np1];
end