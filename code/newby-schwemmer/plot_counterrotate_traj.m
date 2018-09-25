clear

data = load('counterrotate_traj.dat');
x = data(:,2);
y = data(:,3);

figure
plot(x,y, 'k', 'LineWidth', 1)
x = xlabel('x');
y = ylabel('y');
set(x,'Interpreter','latex','fontsize',20)
set(y,'Interpreter','latex','fontsize',20)
axis([-1.5 1.5 -1.5 1.5])
set(gca, 'FontSize', 20)
axis square
% grid on

x = linspace(-1, 1, 1000);
circle_top = sqrt(1-x.^2);
circle_bottom = -sqrt(1-x.^2);

hold on
gray = [0.7 0.7 0.7];
plot(x,circle_top,'color', gray, 'LineWidth', 4)
plot(x,circle_bottom, 'color', gray, 'LineWidth', 4)

x = linspace(-1.5, 1.5, 20);
y = x;
[x, y] = meshgrid(x,y);

omega = 1;
gamma = 15;
c = 4;

u = -omega*y + gamma*x.*(1-x.^2-y.^2) + c*gamma*y.*(x.^2+y.^2-1);
v = omega*x + gamma*y.*(1-x.^2-y.^2) - c*gamma*x.*(x.^2+y.^2-1);

inside = find(x.^2 + y.^2 <=1);
outside = find(x.^2 + y.^2 > 1);

quiver(x(inside), y(inside), u(inside), v(inside), 0.7, 'color', gray, 'LineWidth', 1)
quiver(x(outside), y(outside), u(outside), v(outside), 2, 'color', gray, 'LineWidth', 1)

x = data(:,2);
y = data(:,3);

plot(x,y, 'k', 'LineWidth', 1)
print('-depsc', 'counterrotate-traj')