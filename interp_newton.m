function wx = interp_newton(x, y, x_target)
    w = wyznacz_wspolczynniki_newtona(x, y);
    wx = oblicz_wartosc_wielomianu(x, w, x_target);
end

function w = wyznacz_wspolczynniki_newtona(x, y) 
    if length(x) ~= length(y)
        error('x i y muszą mieć tę samą długość.');
    end
    n = length(x);
    w = [y(1)];
    curr_y = y;

    for i=1:(n-1)
        new_y = [];
        for j=1:(length(curr_y)-1)
            new_y(j) = (curr_y(j+1)-curr_y(j))/(x(j+i)-x(j));
        end
        curr_y = new_y;
        w(i+1) = curr_y(1);
    end
end

function Wx = oblicz_wartosc_wielomianu(x, w, xValue)
    Wx = w(1);
    for i = 2:length(w)
        p = 1;
        for j=2:i
            p = p .* (xValue - x(j-1));
        end
        Wx = Wx + w(i)*p;
    end
end