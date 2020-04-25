function y = stab_fun(x, a)
    y = zeros(size(x));
    for i=1:numel(x)
        y(i) = f(x(i), a(i)) - x(i);
    end
end