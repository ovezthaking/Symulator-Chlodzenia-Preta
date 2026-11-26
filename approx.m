function a = approx(x, y, m, x_target)
    % Skalowanie aby zapobiec ostrzeżeniom o niedokładności
    mu = mean(x);
    sigma = std(x);
    x = (x-mu) / sigma;

    M = zeros(length(x), m+1);
    for i=0:m
        M(:,i+1) = x .^i;
    end
    
    y = y(:);
    a=M \ y;
    if nargin > 3
        x_target = (x_target - mu) / sigma;
        a = oblicz_wielomian(a, x_target);
    end
end

function Fx = oblicz_wielomian(a, x)
    Fx = zeros(size(x));

    for i=0:length(a)-1
        Fx = Fx + a(i+1) .* (x.^i);
    end
end