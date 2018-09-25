function val = g(theta ,r)
sigma = 0.2;
val = r .* (1-r.^2) + r*sigma^2/2.*(cos(theta).^2 - sin(theta).^2);
end