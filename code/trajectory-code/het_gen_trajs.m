function het_gen_trajs(init_cond, seed)

t_end = 100;
dt = 0.001;
t = 0:dt:t_end;

n = 260;

output = NaN(length(t), 2*n+1);
output(:, 1) = t;

data = load('isophase_init_conds.txt');
% data = flipud(data);
x_init = data(init_cond, 1);
y_init = data(init_cond, 2);

output(1, 2:(n+1)) = x_init;
output(1, (n+2):end) = y_init;

rng(seed)

D = 0.01125;

for i = 2:length(t)
    dW_1 = sqrt(dt) * randn(1, n);
    dW_2 = sqrt(dt) * randn(1, n);
    
    output(i, 2:(n+1)) = output(i-1, 2:(n+1)) ...
        + f(output(i-1, 2:(n+1)), output(i-1, (n+2):end))  * dt ...
        + sqrt(2*D) .* dW_1;

    output(i, (n+2):end) = output(i-1, (n+2):end) ...
        + g(output(i-1, 2:(n+1)), output(i-1, (n+2):end))  * dt ...
        + sqrt(2*D) .* dW_2;
    
    % disp(i/length(t)*100)
end

% x = output(:, 22);
% y = output(:, 22+250);
% 
% x = mod(x, 2*pi);
% y = mod(y, 2*pi);
% 
% flip = find(x > 3*pi/2);
% x(flip) = 3*pi/2 - (x(flip) - 3*pi/2);
% 
% flip = find(x < pi/2);
% x(flip) = pi/2 + (pi/2 - x(flip));
% 
% flip = find(y > 3*pi/2);
% y(flip) = 3*pi/2 - (y(flip) - 3*pi/2);
% 
% flip = find(y < pi/2);
% y(flip) = pi/2 + (pi/2 - y(flip));
% 
% x = x - pi;
% y = y - pi;
% 
% x = -x;
% y = -y;

% for i = 1:10:length(x)
%     plot(x(1:i),y(1:i))
%     axis([-pi/2 pi/2 -pi/2 pi/2])
%     axis square
%     pause(0.00001)
% end

save('output', 'output','-v7.3')
end

function val = f(x, y)
alpha = 0.1;
val = cos(x) .* sin(y) + alpha * sin(2*x);
end

function val = g(x, y)
alpha = 0.1;
val = -sin(x) .* cos(y) + alpha * sin(2*y);
end
