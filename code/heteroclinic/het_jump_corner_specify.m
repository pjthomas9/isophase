clear

%% building geometry

% m must be odd!
m = 101; % don't change without checking figure #1
epsilon = 0.1; % don't change without checking figure #1

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

% removing northwest corner
x_all = x_all(2:end);
y_all = y_all(2:end);

z = zeros(length(x_all), 1); % to plot geometry

% geometry plot
figure
plot3(x_all, y_all, z, 'k.', 'MarkerSize', 5)
xlabel('x')
ylabel('y')

% % sorting and plotting boundaries
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

right_jump = find(y_all>=epsilon/2 & y_all<=pi/2 &...
     abs(x_all-h)== min(abs(x_all-h)));
hold on
plot3(x_all(right_jump), y_all(right_jump), z(right_jump),...
    'b.', 'MarkerSize', 10)
set(gca, 'FontSize', 20)
axis square

%% gathering "normal" right and left indices
right_norm = [];
for i = 1:length(y_all)-m
    if y_all(i) == y_all(i+m)
        right_norm = [right_norm i];
    end
end

left_norm = [];
for i = length(y_all):-1:(m+1)
    if y_all(i) == y_all(i-m)
        left_norm = [left_norm i];
    end
end

%% Constructing L-dagger matrix 
D = 0.01125; % noise

N = length(x_all);
L_dagger = NaN(N, N);
for i = 1:N
    temp = zeros(1, N);
    
    temp(i) = -4 * D / h^2; % current spot
    
    up = i - 1; % up index
    down = i + 1; % down index
    
    % find if current spot is on up/down boundary
    check_outer_top = find(outer_top == i); 
    check_inner_bottom = find(inner_bottom == i);
    check_outer_bottom = find(outer_bottom == i);
    check_inner_top = find(inner_top == i);
    
    % inserting "up" coeffecient
    if length(check_outer_top) == 1
        nothing = 0;
    elseif length(check_inner_bottom) == 1
        nothing = 0;
    elseif length(check_outer_bottom) == 1
        temp(up) = 2 * D / h^2;
    elseif length(check_inner_top) == 1
        temp(up) = 2 * D / h^2;
    elseif i >= 2 % northwest corner doesn't exist
        temp(up) = g_het(x_all(i), y_all(i)) / (2*h) + D/h^2;
    end
   
    % inserting "down" coefficient
    if length(check_outer_bottom) == 1
        nothing = 0;
    elseif length(check_inner_top) == 1
        nothing = 0;
    elseif length(check_outer_top) == 1
        temp(down) = 2 * D / h^2;
    elseif length(check_inner_bottom) == 1
        temp(down) = 2 * D / h^2;
    else
        temp(down) = -g_het(x_all(i), y_all(i)) / (2*h) + D/h^2;
    end
    
    % checking if it's in normal right/left index
    check_right_norm = find(right_norm == i);
    check_left_norm = find(left_norm == i);
    
    if length(check_right_norm) == 1
        right = i + m;
    else 
        right = i + (m-r);
    end
    
    if length(check_left_norm) == 1
        left = i - m;
    else 
        left = i - (m-r);
    end
    
    % checking if on left/right boundary
    check_outer_left = find(outer_left == i);
    check_inner_left = find(inner_left == i);
    check_outer_right = find(outer_right == i);
    check_inner_right = find(inner_right == i);  
    
    % filling in coefficients 
    % northwest corner doesn't exist
    if length(check_outer_left) == 0 &&...
            length(check_inner_right) == 0 && i >= (m+1)
        temp(left) = -f_het(x_all(i), y_all(i)) / (2*h) + D/h^2;
    end
    
    if length(check_outer_right) == 0 && length(check_inner_left) == 0
        temp(right) = f_het(x_all(i), y_all(i)) / (2*h) + D/h^2;
    end
    
    if length(check_outer_left) == 1
        temp(right) = 2*D/h^2;
    elseif length(check_outer_right) == 1
        temp(left) = 2*D/h^2;
    elseif length(check_inner_left) == 1
        temp(left) = 2*D/h^2;
    elseif length(check_inner_right) == 1
        temp(right) = 2*D/h^2;
    end

    L_dagger(i, :) = temp;
    disp(i/N*100)
end

%% implementing jump condition in b-vector
b = -ones(N, 1);

T_0 = 100; % specify northwest corner
b(1) = -1 - T_0 * (g_het(x_all(1), y_all(1)) / (2*h) + D/h^2);
b(m) = -1 + T_0 * (f_het(x_all(m), y_all(m)) / (2*h) - D/h^2);

% T_bar = 16.225796508539315; % from simulations
T_bar = 16.233624025276129; % from flux

b(jump) = -1 + T_bar * (f_het(x_all(jump), y_all(jump))/(2*h) + D/h^2);
b(right_jump) = -1 + T_bar * (f_het(x_all(right_jump),...
    y_all(right_jump))/(2*h) - D/h^2);

%% solving system
T = L_dagger \ b;

%% plotting surface
figure
plot3(x_all, y_all, T, 'k.', 'MarkerSize', 5)
hold on
plot3(x_all(outer_top), y_all(outer_top), T(outer_top),...
    'r.', 'MarkerSize', 15)
plot3(x_all(outer_bottom), y_all(outer_bottom), T(outer_bottom),...
    'r.', 'MarkerSize', 15)
plot3(x_all(outer_right), y_all(outer_right), T(outer_right),...
    'r.', 'MarkerSize', 15)
plot3(x_all(outer_left), y_all(outer_left), T(outer_left),...
    'r.', 'MarkerSize', 15)
plot3(x_all(inner_top), y_all(inner_top), T(inner_top),...
    'r.', 'MarkerSize', 15)
plot3(x_all(inner_bottom), y_all(inner_bottom), T(inner_bottom),...
    'r.', 'MarkerSize', 15)
plot3(x_all(inner_right), y_all(inner_right), T(inner_right),...
    'r.', 'MarkerSize', 15)
plot3(x_all(inner_left), y_all(inner_left), T(inner_left),...
    'r.', 'MarkerSize', 15)
for i = 1:5:length(jump)
    plot3([x_all(jump(i)) x_all(right_jump(i))],...
        [y_all(jump(i)) y_all(right_jump(i))],...
        [T(jump(i)) T(right_jump(i))], 'b.-', 'LineWidth', 5)
end

xlabel('x')
ylabel('y')
zlabel('T')
set(gca, 'FontSize', 20)
axis square

%% plotting contourf
T_contour = zeros(m, m);
T_contour(remove) = NaN;

T = [T_0; T];

k = 1;
for i = 1:m^2
    if T_contour(i) == 0
        T_contour(i) = T(k);
        k = k +1;
    end
end

% m = 251; % don't change without checking figure #1
x = linspace(-pi/2, pi/2, m); % column stacked left to right
y = linspace(pi/2, -pi/2, m);

figure
contourf(x, x, flipud(T_contour), 15)
xlabel('x')
ylabel('y')
set(gca, 'FontSize', 20)
axis([-pi/2 pi/2 -pi/2 pi/2])
axis square
colorbar

% save('T_matrix', 'x', 'T_contour')

% C = contourc(x, x, flipud(T_contour), [97.3 97.3]);
% x = C(1, 2:end);
% y = C(2, 2:end);
% 
% iso = [transpose(x) transpose(y)]; 
% [m, n] = size(iso);
% fileID = fopen('isophase_init_conds.txt', 'w');
% 
% hold on
% plot(x(round(m/4.4):end), y(round(m/4.4):m), 'r.-', 'LineWidth', 3)
% 
% num_init_conds = 0;
% for i = round(m/4.4):round(m/35):m
%     fprintf(fileID, '%.11f %.11f \n', iso(i, :));
%     hold on
%     plot(iso(i, 1), iso(i, 2), 'r.', 'MarkerSize', 30)
%     num_init_conds = num_init_conds + 1;
% end
% fclose(fileID);
% 
% 
% fileID = fopen('isophase.txt', 'w');
% for i = round(m/4.4):m
%     fprintf(fileID, '%.11f %.11f \n', iso(i, :));
% end
% fclose(fileID);