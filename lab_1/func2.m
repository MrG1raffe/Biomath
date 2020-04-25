function y = func2(u, v, a)
    y = nthroot(u, 2) .* (1 - a * v.^3);
end