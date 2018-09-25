function val = f(theta, r)
w = 1.99;
k = 1;
sigma = 0.2;
val = w + r*cos(theta) - k*r.^2 + sigma^2/2*cos(theta).*sin(theta);
end