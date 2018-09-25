clear

load('T_matrix.mat');

figure
contourf(X, Y, T, 70)
x = xlabel('x');
y = ylabel('y');
set(x,'Interpreter','latex','fontsize',20)
set(y,'Interpreter','latex','fontsize',20)
set(gca, 'FontSize', 20)
axis([-1.5 1.5 -1.5 1.5])
axis square
colorbar
print('-depsc', 'y-polar-SL-contour')

x = reshape(X, 250^2, 1);
y = reshape(Y, 250^2, 1);
T = reshape(T, 250^2, 1);

figure
plot3(x, y, T, 'k.', 'MarkerSize', 5)

m = 250;

theta = (1:m) * 2 * pi / m;

inner_r = 0.5;
outer_r = 1.5;

r = linspace(outer_r, inner_r, m);

theta_all = repmat(theta, m, 1);
theta_all = reshape(theta_all, m^2, 1);

r_all = diag(r);
r_all = repmat(r_all, m, m);
r_all = diag(r_all);

top = find(r_all == outer_r);
hold on
plot3(x(top), y(top), T(top), 'r.', 'MarkerSize', 20)

bottom = find(r_all == inner_r);
hold on
plot3(x(bottom), y(bottom), T(bottom), 'r.', 'MarkerSize', 20)

left_jump = find(theta_all == theta(1));

right_jump = find(theta_all == theta(end));

for i = 1:30:length(left_jump)
    plot3([x(left_jump(i)) x(right_jump(i))], [y(left_jump(i)) y(right_jump(i))], [T(left_jump(i)) T(right_jump(i))], 'b.-', 'LineWidth', 5)
end

x = xlabel('x');
y = ylabel('y');
z = zlabel('T');
set(x,'Interpreter','latex','fontsize',20)
set(y,'Interpreter','latex','fontsize',20)
set(z,'Interpreter','latex','fontsize',20)
set(gca, 'FontSize', 20)
axis square
print('-depsc', 'y-polar-SL-surface')