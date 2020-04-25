function mu = lyap(x0, a)
    n = 1000;
    mu = zeros(n, 1);
    mu(1) = x0;
    for i = 2:n
        mu(i) = f(mu(i-1), a);
    end
    for i = 1:n
        mu(i) = f_der(mu(i), a);
    end
    mu = sum(log(mu)) / n;
end