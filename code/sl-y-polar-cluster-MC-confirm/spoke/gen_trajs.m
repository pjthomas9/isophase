function gen_trajs(init_cond, seed)

t_end = 100;
dt = 0.001;
t = 0:dt:t_end;

n = 300;

output = NaN(length(t), 2*n+1);
output(:, 1) = t;

data = load('spoke_init_conds.txt');
% data = flipud(data);
theta_init = data(init_cond, 1);
r_init = data(init_cond, 2);

output(1, 2:(n+1)) = theta_init;
output(1, (n+2):end) = r_init;

rng(seed)

for i = 2:length(t)
    dW = sqrt(dt) * randn(1, n);
    
    output(i, 2:(n+1)) = output(i-1, 2:(n+1)) ...
        + f(output(i-1, 2:(n+1)), output(i-1, (n+2):end))  * dt ...
        + B1(output(i-1, 2:(n+1)), output(i-1, (n+2):end)) .* dW;

    output(i, (n+2):end) = output(i-1, (n+2):end) ...
        + g(output(i-1, 2:(n+1)), output(i-1, (n+2):end))  * dt ...
        + B2(output(i-1, 2:(n+1)), output(i-1, (n+2):end)) .* dW;
    
    disp(i/length(t)*100)
end

% theta = output(:, 2);
% r = output(:, 2+n);
% 
% x = r .* cos(theta);
% y = r .* sin(theta);
% 
% for i = 1:10:length(x)
%     plot(x(1:i),y(1:i))
%     axis([-1.5 1.5 -1.5 1.5])
%     axis square
%     pause(0.00001)
% end

save('output', 'output','-v7.3')
end

function val = f(theta, r)
w=1.99;
k=1;
sigma=0.2;
val = w+r.*cos(theta)-k*r.^2+0.5*(sigma^2*cos(theta).*sin(theta));
end

function val = g(theta, r)
w=1.99;
k=1;
sigma=0.2;
val = r.*(1-r.^2)+0.5*sigma^2*r.*(cos(theta).^2 - sin(theta).^2);
end

function val = B1(theta, r)
w=1.99;
k=1;
sigma=0.2;
val = sigma*sin(theta);
end

function val = B2(theta, r)
w=1.99;
k=1;
sigma=0.2;
val = sigma*r.*cos(theta);
end
