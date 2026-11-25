function a = approx(x, y, m, x_target)
    M = zeros(length(x), m+1);
    for i=0:m
        M(:,i+1) = x .^i;
    end
    
    y = y(:);
    a=(M' * M) \ (M' * y);
    if nargin > 3
        a = oblicz_wielomian(a, x_target);
    end
end

function Fx = oblicz_wielomian(a, x)
    Fx = zeros(size(x));

    for i=0:length(a)-1
        Fx = Fx + a(i+1) .* (x.^i);
    end
end