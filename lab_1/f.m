function y = f(x, a)
    y = nthroot(x, 2) .* (1 - a * x.^3);
end