function val = BBT12(theta, r)
sigma = 0.2;
val = sigma^2 * r .* sin(theta) .* cos(theta);
end