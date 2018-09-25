function calc_frts(job_ID)

data = load('spoke_init_conds.txt');
slope = data(1,2) / data(1,1);

x_spoke = linspace(data(1,1), 0, 100);
y_spoke = slope * x_spoke;

load('output.mat');
t = output(:,1);

n = 260;
first_return_times = [];

for i = 1:n
    x = output(:, i+1);
    y = output(:, i+1+n);

    x = mod(x, 2*pi);
    y = mod(y, 2*pi);

    flip = find(x > 3*pi/2);
    x(flip) = 3*pi/2 - (x(flip) - 3*pi/2);

    flip = find(x < pi/2);
    x(flip) = pi/2 + (pi/2 - x(flip));

    flip = find(y > 3*pi/2);
    y(flip) = 3*pi/2 - (y(flip) - 3*pi/2);

    flip = find(y < pi/2);
    y(flip) = pi/2 + (pi/2 - y(flip));

    x = x - pi;
    y = y - pi;

    x = -x;
    y = -y;
    
    theta = atan(y./x);
    
    d_theta = [0; diff(theta)];
    
    check = find(d_theta > 2);
    d_theta(check) = d_theta(check) - pi;
    
    check = find(d_theta < -2);
    d_theta(check) = d_theta(check) + pi;
    
    theta_sum = cumsum(d_theta);
    
    stop = find(theta_sum <= -2*pi);
    if length(stop) >= 1
        stop = stop(1);
        first_return_times = [first_return_times t(stop)];
        
        
        figure(1)
        clf
        subplot(1,2,1)
        plot(x_spoke, y_spoke, 'r', 'LineWidth', 2)
        hold on
        plot(x(1:stop), y(1:stop), 'b')
        axis([-pi/2 pi/2 -pi/2 pi/2])
        axis square
        title(num2str(i))
        
        subplot(1,2,2)
        plot(x_spoke, y_spoke, 'r', 'LineWidth', 2)
        hold on
        plot(x(1:stop), y(1:stop), 'b')
        axis([x(stop)-0.1 x(stop)+0.1 y(stop)-0.1 y(stop)+0.1])
        axis square
        title(num2str(t(stop)))
   
        pause(0.5)
        
    end
end

save(['frts_' num2str(job_ID)], 'first_return_times')
end
