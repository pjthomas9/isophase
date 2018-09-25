function calc_frts(job_ID)

spoke_data = load('spoke_init_conds.txt');
x_spoke = spoke_data(:,2) .* cos(spoke_data(:,1));
y_spoke = spoke_data(:,2) .* sin(spoke_data(:,1));


load('output.mat');
t = output(:,1);

n = 300;
first_return_times = [];

for i = 1:n
    theta = output(:, i+1);
    r = output(:, i+1+n);
    
    x = r.*cos(theta);
    y = r.*sin(theta);
    
    d_theta = [0; diff(theta)];
    
    theta_sum = cumsum(d_theta);
    
    stop = find(theta_sum >= 2*pi);
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
   
        pause(0.1)
        
    end
end
mean(first_return_times)

save(['frts_' num2str(job_ID)], 'first_return_times')
end
