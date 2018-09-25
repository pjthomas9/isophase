clear

m = 50;

x = linspace(-pi/2, pi/2, m); 
y = linspace(pi/2, -pi/2, m);

x_all = repmat(x, m, 1); 
x_all = reshape(x_all, m^2, 1);

y_all = diag(y);
y_all = repmat(y_all, m, m);
y_all = diag(y_all);

T = atan(y_all ./ x_all);
shift = find(x_all > 0);
T(shift) = T(shift) + pi;

figure
plot3(x_all, y_all, T, 'k.')
x = xlabel('x');
y = ylabel('y');
z = zlabel('z');
set(x,'Interpreter','latex','fontsize',20)
set(y,'Interpreter','latex','fontsize',20)
set(z,'Interpreter','latex','fontsize',20)

hold on
plot3(x_all, y_all, T+2*pi, 'k.')
plot3(x_all, y_all, T+4*pi, 'k.')

set(gca, 'FontSize', 20)
axis square

print('-depsc', 'time-crystal')