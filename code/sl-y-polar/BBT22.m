function val = BBT22(theta, r)
sigma = 0.2;
val = sigma^2 * r.^2 .* cos(theta).^2;
end