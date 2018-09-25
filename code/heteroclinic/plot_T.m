clear

load('T_matrix.mat');

figure
contourf(x, x, flipud(T_contour), 20)
x = xlabel('x');
y = ylabel('y');
set(x,'Interpreter','latex','fontsize',20)
set(y,'Interpreter','latex','fontsize',20)
set(gca, 'FontSize', 20)
axis([-pi/2 pi/2 -pi/2 pi/2])
axis square
colorbar
print('-depsc', 'het-T-contour')

% m must be odd!
m = 251; % don't change without checking figure #1
epsilon = 0.05; % don't change without checking figure #1

x = linspace(-pi/2, pi/2, m); % column stacked left to right
y = linspace(pi/2, -pi/2, m);

h = diff(x); % discretization size
h = h(1);

x_all = repmat(x, m, 1); % forming vector of ALL x coordinates in order
x_all = reshape(x_all, m^2, 1);

y_all = diag(y); % forming vector of ALL y coordinates in order
y_all = repmat(y_all, m, m);
y_all = diag(y_all);

remove = find(x_all>-epsilon/2 & x_all<epsilon/2 & y_all>-epsilon/2 ...
    & y_all< epsilon/2);
r = sqrt(length(remove));
x_all(remove) = [];
y_all(remove) = [];

outer_top = find(y_all == pi/2);

outer_bottom = find(y_all == -pi/2);

outer_right = find(x_all == pi/2);

outer_left = find(x_all == -pi/2);

inner_top = find(x_all>=-epsilon/2 & x_all<=epsilon/2 &...
    abs(y_all-epsilon/2) <= h);

inner_bottom = find(x_all>=-epsilon/2 & x_all<=epsilon/2 &...
    abs(y_all+epsilon/2) <= h);

inner_right = find(y_all>=-epsilon/2 & y_all<=epsilon/2 &...
    abs(x_all-epsilon/2) <= h);

inner_left = find(y_all>=-epsilon/2 & y_all<=epsilon/2 &...
    abs(x_all+epsilon/2) <= h);

jump = find(y_all>=epsilon/2 & y_all<=pi/2 & x_all == 0);

right_jump = find(y_all>=epsilon/2 & y_all<=pi/2 &...
    abs(x_all-h)== min(abs(x_all-h)));

T = reshape(T_contour, m^2, 1);
remove = find(T > 50);
T = T(remove);

figure
plot3(x_all, y_all, T, 'k.', 'MarkerSize', 0.05)
hold on
plot3(x_all(outer_top), y_all(outer_top), T(outer_top),...
    'r.', 'MarkerSize', 12)
plot3(x_all(outer_bottom), y_all(outer_bottom), T(outer_bottom),...
    'r.', 'MarkerSize', 12)
plot3(x_all(outer_right), y_all(outer_right), T(outer_right),...
    'r.', 'MarkerSize', 12)
plot3(x_all(outer_left), y_all(outer_left), T(outer_left),...
    'r.', 'MarkerSize', 12)
plot3(x_all(inner_top), y_all(inner_top), T(inner_top),...
    'r.', 'MarkerSize', 12)
plot3(x_all(inner_bottom), y_all(inner_bottom), T(inner_bottom),...
    'r.', 'MarkerSize', 12)
plot3(x_all(inner_right), y_all(inner_right), T(inner_right),...
    'r.', 'MarkerSize', 12)
plot3(x_all(inner_left), y_all(inner_left), T(inner_left),...
    'r.', 'MarkerSize', 12)
for i = 1:5:length(jump)
    plot3([x_all(jump(i)) x_all(right_jump(i))],...
        [y_all(jump(i)) y_all(right_jump(i))],...
        [T(jump(i)) T(right_jump(i))], 'b.-', 'LineWidth', 2)
end

x = xlabel('x');
y = ylabel('y');
z = zlabel('T');
set(x,'Interpreter','latex','fontsize',20)
set(y,'Interpreter','latex','fontsize',20)
set(z,'Interpreter','latex','fontsize',20)
set(gca, 'FontSize', 20)
axis square
print('-depsc', 'het-T-surface')