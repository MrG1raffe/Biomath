function y = f_der(u, a)
    y = (0.5 * (1 - a * u.^3) ./ nthroot(u, 2) - 3 * a * u .^ 2.5); 
end