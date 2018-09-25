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

print('-depsc', 'het-traj')




