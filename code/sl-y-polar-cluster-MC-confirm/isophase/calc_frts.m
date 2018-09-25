function calc_frts(job_ID)

isophase = load('isophase.txt');
% isophase = flipud(isophase);

iso_x = isophase(:, 1);
iso_y = isophase(:, 2);

load('output.mat');
t = output(:,1);

n = 300;
first_return_times = [];

for i = 1:n
    theta = output(:, i+1);
    r = output(:, i+1+n);

    x = r .* cos(theta);
    y = r .* sin(theta);

    test = NaN(length(x), 1);
    dist = NaN(length(iso_x), 1);

    for j = 1:length(test)
        dist = sqrt((x(j) - iso_x).^2 + (y(j) - iso_y).^2);
        test(j) = min(dist);
    end
    
    d_theta = [0; diff(theta)];
    
    theta_sum = cumsum(d_theta);
    not_far_enough = find(theta_sum <= 5);
    
    too_far = find(theta_sum >= 2*pi*2);
    
    test(not_far_enough) = 5;
    test(too_far) = 5;
    
    stop = find(test < 0.01);
    if length(stop) >= 1
        stop = stop(1);
        first_return_times = [first_return_times t(stop)];
        
        
%         figure(1)
%         clf
%         subplot(1,2,1)
%         plot(iso_x, iso_y, 'r', 'LineWidth', 2)
%         hold on
%         plot(x(1:stop), y(1:stop), 'b')
%         
%         circ_theta = 0:0.01:2*pi;
%         plot(0.5*cos(circ_theta), 0.5*sin(circ_theta), 'r', 'LineWidth', 3)
%         
%         axis([-1.5 1.5 -1.5 1.5])
%         axis square
%         title(num2str(i))
%         
%         subplot(1,2,2)
%         plot(iso_x, iso_y, 'r', 'LineWidth', 2)
%         hold on
%         plot(x(1:stop), y(1:stop), 'b')
%         axis([x(stop)-0.1 x(stop)+0.1 y(stop)-0.1 y(stop)+0.1])
%         axis square
%         title(num2str(t(stop)))
%    
%         pause(0.01)
        
    end
end

% format long
% mean(first_return_times)

save(['frts_' num2str(job_ID)], 'first_return_times')
end
