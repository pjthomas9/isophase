function val = f_het(x, y)
alpha = 0.1;
val = cos(x) .* sin(y) + alpha * sin(2*x);
end