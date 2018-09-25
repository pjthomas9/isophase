function val = g_het(x, y)
alpha = 0.1;
val = -sin(x) .* cos(y) + alpha * sin(2*y);
end