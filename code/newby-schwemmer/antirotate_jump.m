clear

n = 80;

inner_r = 0.5;
outer_r = 1.5;

r = linspace(inner_r, outer_r, n);
dr = diff(r);
dr = dr(1);

theta = (1:n) * 2 * pi / n;
d_theta = diff(theta);
d_theta = d_theta(1);

r_all = repmat(r, n, 1);
r_all = reshape(r_all, n^2, 1);

theta_all = diag(theta);
theta_all = repmat(theta_all, n, n);
theta_all = diag(theta_all);

z = zeros(length(r_all), 1);
plot3(r_all.*cos(theta_all), r_all.*sin(theta_all), z, 'k.', 'MarkerSize', 1)

outer_ring = find(r_all == outer_r);
hold on
plot3(r_all(outer_ring).*cos(theta_all(outer_ring)), r_all(outer_ring).*sin(theta_all(outer_ring)), z(outer_ring), 'r.', 'MarkerSize', 10)

inner_ring = find(r_all == inner_r);
hold on
plot3(r_all(inner_ring).*cos(theta_all(inner_ring)), r_all(inner_ring).*sin(theta_all(inner_ring)), z(inner_ring), 'r.', 'MarkerSize', 5)

left_jump = find(theta_all == theta(1));
hold on
plot3(r_all(left_jump).*cos(theta_all(left_jump)), r_all(left_jump).*sin(theta_all(left_jump)), z(left_jump), 'b.', 'MarkerSize', 10)

right_jump = find(theta_all == theta(end));
hold on
plot3(r_all(right_jump).*cos(theta_all(right_jump)), r_all(right_jump).*sin(theta_all(right_jump)), z(right_jump), 'b.', 'MarkerSize', 10)
xlabel('x')
ylabel('y')
axis square
% set(gca, 'FontSize', 20)

D = 0.19845;
gamma = 15;
w = 1;
c = -15;

N = length(r_all);
L_dagger = NaN(N, N);

for i = 1:N
    temp = zeros(1, N);
    
    temp(i) = D * (r_all(i)^(-2) * (-2 / d_theta^2) - 2 / dr^2);
    
    in = i - n;
    out = i + n;
    
    if mod(i, n) >= 2
        clockwise = i - 1;
        counter_clockwise = i + 1;
    elseif mod(i, n) == 1
        clockwise = i - 1 + n;
        counter_clockwise = i + 1;
    elseif mod(i, n) == 0
        clockwise = i - 1;
        counter_clockwise = i + 1 - n;
    end
    
    check_inner_ring = find(inner_ring == i);
    check_outer_ring = find(outer_ring == i);
    
    if length(check_inner_ring) == 1
        temp(out) = D * (2 / dr^2);
    elseif length(check_outer_ring) == 1
        temp(in) = D * (2 / dr^2);
    else
        temp(in) = (-gamma * r_all(i) * (1-r_all(i)^2) - ...
            D* r_all(i)^(-1)) / (2 * dr) + D / dr^2;
        temp(out) = (gamma * r_all(i) * (1-r_all(i)^2) + ...
            D * r_all(i)^(-1)) / (2 * dr) + D / dr^2;
    end

	temp(clockwise) = -(w +c*gamma*w*(1-r_all(i))^2) / (2 * d_theta) + D / (r_all(i) * d_theta)^2;
	temp(counter_clockwise) = (w +c*gamma*w*(1-r_all(i))^2) / (2 * d_theta) + D / (r_all(i) * d_theta)^2;
    
    L_dagger(i, :) = temp;
    disp(i/N*100)
end

T_bar = -11.890000934672861;

b = -ones(N, 1);
b(left_jump) = -1 + T_bar * ((w +c*gamma*w*(1-r_all(left_jump)).^2) / (2*d_theta) - D ./ (r_all(left_jump).^2 * d_theta^2));
b(right_jump) = -1 + T_bar * ((w +c*gamma*w*(1-r_all(right_jump)).^2) / (2*d_theta) + D ./ (r_all(right_jump).^2 * d_theta^2));

T = L_dagger \ b;
T = T - min(T);

figure
plot3(r_all.*cos(theta_all), r_all.*sin(theta_all), T, 'k.', 'MarkerSize', 1)
hold on
plot3(r_all(outer_ring).*cos(theta_all(outer_ring)), r_all(outer_ring).*sin(theta_all(outer_ring)), T(outer_ring), 'r.', 'MarkerSize', 15)
plot3(r_all(inner_ring).*cos(theta_all(inner_ring)), r_all(inner_ring).*sin(theta_all(inner_ring)), T(inner_ring), 'r.', 'MarkerSize', 15)

x_all = r_all .* cos(theta_all);
y_all = r_all .* sin(theta_all);

for i = 1:15:length(left_jump)
    plot3([x_all(left_jump(i)) x_all(right_jump(i))], [y_all(left_jump(i)) y_all(right_jump(i))], [T(left_jump(i)) T(right_jump(i))], 'b.-', 'LineWidth', 3)
end
x = xlabel('x');
y = ylabel('y');
z = zlabel('T');
set(x,'Interpreter','latex','fontsize',20)
set(y,'Interpreter','latex','fontsize',20)
set(z,'Interpreter','latex','fontsize',20)
axis square
set(gca, 'FontSize', 20)
print('-depsc', 'antirotate-surface')

X = r_all .* cos(theta_all);
X = reshape(X, n, n);

Y = r_all .* sin(theta_all);
Y = reshape(Y, n, n);

T = reshape(T, n, n);
figure
contourf(X, Y, T, 30)
x = xlabel('x');
y = ylabel('y');
set(x,'Interpreter','latex','fontsize',20)
set(y,'Interpreter','latex','fontsize',20)
axis square
colorbar
set(gca, 'FontSize', 20)
print('-depsc', 'antirotate-contour')