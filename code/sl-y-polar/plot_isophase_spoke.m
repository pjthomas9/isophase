clear

data = load('y-polar-SL-ex-traj.dat');
x = data(:,4);
y = data(:,5);

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
x = data(:,2) .* cos(data(:,1));
y = data(:,2) .* sin(data(:,1));
plot(x, y, 'b.', 'MarkerSize', 40)

theta = data(10, 1);

r = data(:,2);

x = r * cos(theta);
y = r * sin(theta);

plot(x,y, 'r.', 'MarkerSize', 30)

r = [r; 0];
x = r * cos(theta);
y = r * sin(theta);
plot(x,y, 'r-', 'LineWidth', 2)

isophase = load('isophase.txt');
plot(isophase(:,1), isophase(:,2), 'b', 'LineWidth', 4)

print('-depsc', 'y-polar-SL-spoke-isophase')


