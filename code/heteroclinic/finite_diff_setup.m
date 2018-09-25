clear

%% building geometry

% m must be odd!
m = 101; % don't change without checking figure #1
epsilon = 0.3; % don't change without checking figure #1

x = linspace(-pi/2, pi/2, m); % column stacked left to right
y = linspace(pi/2, -pi/2, m);

h = diff(x); % discretization size
h = h(1);

x_all = repmat(x, m, 1); % forming vector of ALL x coordinates in order
x_all = reshape(x_all, m^2, 1);

y_all = diag(y); % forming vector of ALL y coordinates in order
y_all = repmat(y_all, m, m);
y_all = diag(y_all);

% removing center square
remove = find(x_all>-epsilon/2 & x_all<epsilon/2 & y_all>-epsilon/2 ...
    & y_all< epsilon/2);
r = sqrt(length(remove));
x_all(remove) = [];
y_all(remove) = [];

z = zeros(length(x_all), 1); % to plot geometry

% geometry plot
figure
plot3(x_all, y_all, z, 'k.', 'MarkerSize', 5)
x = xlabel('x');
y = ylabel('y');
set(x,'Interpreter','latex','fontsize',20)
set(y,'Interpreter','latex','fontsize',20)

% sorting and plotting boundaries
outer_top = find(y_all == pi/2);
hold on
plot3(x_all(outer_top), y_all(outer_top), z(outer_top), 'r.',...
    'MarkerSize', 10)

outer_bottom = find(y_all == -pi/2);
hold on
plot3(x_all(outer_bottom), y_all(outer_bottom), z(outer_bottom),...
    'r.', 'MarkerSize', 10)

outer_right = find(x_all == pi/2);
hold on
plot3(x_all(outer_right), y_all(outer_right), z(outer_right),...
    'r.', 'MarkerSize', 10)

outer_left = find(x_all == -pi/2);
hold on
plot3(x_all(outer_left), y_all(outer_left), z(outer_left), 'r.',...
    'MarkerSize', 10)

inner_top = find(x_all>=-epsilon/2 & x_all<=epsilon/2 &...
    abs(y_all-epsilon/2) <= h);
hold on
plot3(x_all(inner_top), y_all(inner_top), z(inner_top), 'r.',...
    'MarkerSize', 10)

inner_bottom = find(x_all>=-epsilon/2 & x_all<=epsilon/2 &...
    abs(y_all+epsilon/2) <= h);
hold on
plot3(x_all(inner_bottom), y_all(inner_bottom), z(inner_bottom),...
    'r.', 'MarkerSize', 10)

inner_right = find(y_all>=-epsilon/2 & y_all<=epsilon/2 &...
    abs(x_all-epsilon/2) <= h);
hold on
plot3(x_all(inner_right), y_all(inner_right), z(inner_right),...
    'r.', 'MarkerSize', 10)

inner_left = find(y_all>=-epsilon/2 & y_all<=epsilon/2 &...
    abs(x_all+epsilon/2) <= h);
hold on
plot3(x_all(inner_left), y_all(inner_left), z(inner_left),...
    'r.', 'MarkerSize', 10)

jump = find(y_all>=epsilon/2 & y_all<=pi/2 & x_all == 0);
hold on
plot3(x_all(jump), y_all(jump), z(jump), 'b.', 'MarkerSize', 10)


set(gca, 'FontSize', 20)
axis square
view([0,90])

print('-depsc', 'het-finite-diff')