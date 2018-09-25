clear

data = load('het_traj.dat');
x = data(:,2);
y = data(:,3);

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

figure
plot(x,y, 'k')
x = xlabel('x');
y = ylabel('y');
set(x,'Interpreter','latex','fontsize',20)
set(y,'Interpreter','latex','fontsize',20)
axis([-pi/2 pi/2 -pi/2 pi/2])
set(gca, 'FontSize', 20)
axis square
grid on

data = load('isophase_init_conds.txt');
hold on
plot(data(:,1), data(:,2), 'b.', 'MarkerSize', 40)
middle = data(10, :);

theta = atan(middle(2) / middle(1)) + pi;

r = sqrt(data(:,1).^2 + data(:,2).^2);

x = r * cos(theta);
y = r * sin(theta);

plot(x,y, 'r.', 'MarkerSize', 40)

r = [r; 0];
x = r * cos(theta);
y = r * sin(theta);
plot(x,y, 'r-', 'LineWidth', 4)

isophase = load('isophase.txt');
plot(isophase(:,1), isophase(:,2), 'b', 'LineWidth', 4)

print('-depsc', 'het-spoke-isophase')


