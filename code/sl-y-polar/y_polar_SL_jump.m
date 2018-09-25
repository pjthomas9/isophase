clear

m = 100;

theta = (1:m) * 2 * pi / m;
d_theta = diff(theta);
d_theta = d_theta(1);

inner_r = 0.5;
outer_r = 1.5 ;

r = linspace(outer_r, inner_r, m);
dr = diff(r);
dr = abs(dr(1));

theta_all = repmat(theta, m, 1);
theta_all = reshape(theta_all, m^2, 1);

r_all = diag(r);
r_all = repmat(r_all, m, m);
r_all = diag(r_all);

z = zeros(size(r_all));

figure
plot3(theta_all, r_all, z, 'k.', 'MarkerSize', 5)
xlabel('\theta')
ylabel('r')

top = find(r_all == outer_r);
hold on
plot3(theta_all(top), r_all(top), z(top), 'r.', 'MarkerSize', 10)

bottom = find(r_all == inner_r);
hold on
plot3(theta_all(bottom), r_all(bottom), z(bottom), 'r.', 'MarkerSize', 10)

left_jump = find(theta_all == theta(1));
hold on
plot3(theta_all(left_jump), r_all(left_jump), z(left_jump), 'b.', 'MarkerSize', 10)

right_jump = find(theta_all == theta(end));
hold on
plot3(theta_all(right_jump), r_all(right_jump), z(right_jump), 'b.', 'MarkerSize', 10)
set(gca, 'FontSize', 20)
axis square
axis([0 2*pi 0.5 1.5 -1 1])

w = 1.99;
k = 1;
sigma = 0.2;

N = length(r_all);
extra = length(top);
L_dagger = NaN(N+2*extra, N);

b = -ones(N+2*extra, 1);
b(N+1:end) = 0;
% T_bar = 20.739357674673091; % from sims
T_bar = 20.929154018932547; % from flux

counter = 1;

for i = 1:N
    temp = zeros(1, N);
    temp_bc = temp;

    if r_all(i) ~= inner_r && r_all(i) ~= outer_r
        temp(i) = -BBT11(theta_all(i), r_all(i))/d_theta^2 - BBT22(theta_all(i), r_all(i))/dr^2;

        up = i - 1;
        down = i + 1;
        
        temp(up) = g(theta_all(i) , r_all(i)) / (2*dr) + 0.5*BBT22(theta_all(i), r_all(i)) / dr^2;
        temp(down) = -g(theta_all(i) , r_all(i)) / (2*dr) + 0.5*BBT22(theta_all(i), r_all(i)) / dr^2;

        left = i - m;
        right = i + m;
    
        if theta_all(i) == theta(1)
            left = i + m*(m-1);
        elseif theta_all(i) == theta(end)
            right = i - m*(m-1);
        end

        temp(right) = f(theta_all(i), r_all(i)) / (2*d_theta) + 0.5*BBT11(theta_all(i), r_all(i)) / d_theta^2;
        temp(left) = -f(theta_all(i), r_all(i)) / (2*d_theta) + 0.5*BBT11(theta_all(i), r_all(i)) / d_theta^2;

        if theta_all(i) == theta(1)
            b(i) = b(i) + T_bar * (f(theta_all(i), r_all(i)) / (2*d_theta) - 0.5*BBT11(theta_all(i), r_all(i)) / d_theta^2);
        end

        if theta_all(i) == theta(end)
            b(i) = b(i) + T_bar * (f(theta_all(i), r_all(i)) / (2*d_theta) + 0.5*BBT11(theta_all(i), r_all(i)) / d_theta^2);
        end
    
        NW = i - m - 1;
        SW = i - m + 1;
        NE = i + m - 1;
        SE = i + m + 1;
    
        if theta_all(i) == theta(1)
            NW = i + m*(m-1) - 1; 
            SW = i + m*(m-1) + 1;
        elseif theta_all(i) == theta(end)
            NE = i - m*(m-1) - 1;
            SE = i - m*(m-1) + 1;
        end

        temp(NE) = BBT12(theta_all(i), r_all(i)) / (4*d_theta*dr);
        temp(SE) = -BBT12(theta_all(i), r_all(i)) / (4*d_theta*dr);
        temp(NW) = -BBT12(theta_all(i), r_all(i)) / (4*d_theta*dr);
        temp(SW) = BBT12(theta_all(i), r_all(i)) / (4*d_theta*dr);

        L_dagger(i, :) = temp;
        disp(i/N*100)

    elseif r_all(i) == outer_r
        temp(i) = g(theta_all(i), r_all(i)) * (1/dr) * (3/2) ...
                    - BBT11(theta_all(i), r_all(i)) / d_theta^2 ...
                    + 0.5 * BBT22(theta_all(i), r_all(i)) / dr^2;

        down = i + 1;
        down2 = i + 2;
        
        temp(down) = g(theta_all(i), r_all(i)) * (1/dr) * (-2) ...
                        - BBT22(theta_all(i), r_all(i)) / dr^2;

        temp(down2) = g(theta_all(i), r_all(i)) * (1/dr) * (1/2) ...
                        + 0.5 * BBT22(theta_all(i), r_all(i)) / dr^2;


        left = i - m;
        right = i + m;
    
        if theta_all(i) == theta(1)
            left = i + m*(m-1);
        elseif theta_all(i) == theta(end)
            right = i - m*(m-1);
        end

        temp(right) = f(theta_all(i), r_all(i)) / (2*d_theta) ...
                        + 0.5 * BBT11(theta_all(i), r_all(i)) / d_theta^2 ...
                        + (1/dr) * BBT12(theta_all(i), r_all(i)) * (3/2) / (2 * d_theta);
        temp(left) = -f(theta_all(i), r_all(i)) / (2*d_theta) ...
                        + 0.5 * BBT11(theta_all(i), r_all(i)) / d_theta^2 ...
                        - (1/dr) * BBT12(theta_all(i), r_all(i)) * (3/2) / (2 * d_theta);

        if theta_all(i) == theta(1)
            b(i) = b(i) + T_bar * (f(theta_all(i), r_all(i)) / (2*d_theta)...
                            - 0.5 * BBT11(theta_all(i), r_all(i)) / d_theta^2 ...
                            + (1/dr) * BBT12(theta_all(i), r_all(i)) * (3/2) / (2 * d_theta));
        end

        if theta_all(i) == theta(end)
            b(i) = b(i) + T_bar * (f(theta_all(i), r_all(i)) / (2*d_theta)...
                            + 0.5 * BBT11(theta_all(i), r_all(i)) / d_theta^2 ...
                            + (1/dr) * BBT12(theta_all(i), r_all(i)) * (3/2) / (2 * d_theta));
        end
                    
        
        DL = i - m + 1;
        DR = i + m + 1;
        DDL = i - m + 2;
        DDR = i + m + 2;
    
        if theta_all(i) == theta(1)
            DL = i + m*(m-1) + 1;
            DDL = DL + 1;
        elseif theta_all(i) == theta(end)
            DR = i - m*(m-1) + 1;
            DDR = DR + 1;
        end

        temp(DR) = -BBT12(theta_all(i), r_all(i)) / dr / d_theta;
        temp(DL) = BBT12(theta_all(i), r_all(i)) / dr / d_theta;

        temp(DDR) = BBT12(theta_all(i), r_all(i)) / dr / 4 / d_theta;
        temp(DDL) = -BBT12(theta_all(i), r_all(i)) / dr / 4 / d_theta;


        if theta_all(i) == theta(1)
            b(i) = b(i) + T_bar * (-BBT12(theta_all(i), r_all(i)) / dr / d_theta...
                                       +BBT12(theta_all(i), r_all(i)) / dr / 4 / d_theta);
        end

        if theta_all(i) == theta(end)
            b(i) = b(i) + T_bar * (-BBT12(theta_all(i), r_all(i)) / dr / d_theta...
                                       +BBT12(theta_all(i), r_all(i)) / dr / 4 / d_theta);
        end

        temp_bc(i) = BBT22(theta_all(i), r_all(i)) / dr * (3/2);
        temp_bc(right) = BBT12(theta_all(i), r_all(i)) / 2 / d_theta;
        temp_bc(left) = -BBT12(theta_all(i), r_all(i)) / 2 / d_theta;
        temp_bc(down) = BBT22(theta_all(i), r_all(i)) / dr * (-2);
        temp_bc(down2) = BBT22(theta_all(i), r_all(i)) / dr * (1/2);

        if theta_all(i) == theta(1)
            b(N+counter) = b(N+counter) + T_bar * (BBT12(theta_all(i), r_all(i)) / 2 / d_theta);
        end

        if theta_all(i) == theta(end)
            b(N+counter) = b(N+counter) + T_bar * (BBT12(theta_all(i), r_all(i)) / 2 / d_theta);
        end

        L_dagger(i, :) = temp;
        L_dagger(N+counter, :) = temp_bc;
        counter = counter + 1;
        disp(i/N*100)
    elseif r_all(i) == inner_r
        temp(i) = g(theta_all(i), r_all(i)) * (1/dr) * (-3/2) ...
                    - BBT11(theta_all(i), r_all(i)) / d_theta^2 ...
                    + 0.5 * BBT22(theta_all(i), r_all(i)) / dr^2;

        up = i - 1;
        up2 = i - 2;
        
        temp(up) = g(theta_all(i), r_all(i)) * (1/dr) * (2) ...
                        - BBT22(theta_all(i), r_all(i)) / dr^2;

        temp(up2) = g(theta_all(i), r_all(i)) * (1/dr) * (-1/2) ...
                        + 0.5 * BBT22(theta_all(i), r_all(i)) / dr^2;


        left = i - m;
        right = i + m;
    
        if theta_all(i) == theta(1)
            left = i + m*(m-1);
        elseif theta_all(i) == theta(end)
            right = i - m*(m-1);
        end

        temp(right) = f(theta_all(i), r_all(i)) / (2*d_theta) ...
                        + 0.5 * BBT11(theta_all(i), r_all(i)) / d_theta^2 ...
                        + (1/dr) * BBT12(theta_all(i), r_all(i)) * (-3/2) / (2 * d_theta);
        temp(left) = -f(theta_all(i), r_all(i)) / (2*d_theta) ...
                        + 0.5 * BBT11(theta_all(i), r_all(i)) / d_theta^2 ...
                        - (1/dr) * BBT12(theta_all(i), r_all(i)) * (-3/2) / (2 * d_theta);

        if theta_all(i) == theta(1)
            b(i) = b(i) + T_bar * (f(theta_all(i), r_all(i)) / (2*d_theta)...
                            - 0.5 * BBT11(theta_all(i), r_all(i)) / d_theta^2 ...
                            - (1/dr) * BBT12(theta_all(i), r_all(i)) * (3/2) / (2 * d_theta));
        end

        if theta_all(i) == theta(end)
            b(i) = b(i) + T_bar * (f(theta_all(i), r_all(i)) / (2*d_theta)...
                            + 0.5 * BBT11(theta_all(i), r_all(i)) / d_theta^2 ...
                            - (1/dr) * BBT12(theta_all(i), r_all(i)) * (3/2) / (2 * d_theta));
        end
                    
        
        UL = i - m - 1;
        UR = i + m - 1;
        UUL = i - m - 2;
        UUR = i + m - 2;
    
        if theta_all(i) == theta(1)
            UL = i + m*(m-1) - 1;
            UUL = UL - 1;
        elseif theta_all(i) == theta(end)
            UR = i - m*(m-1) - 1;
            UUR = DR - 1;
        end

        temp(UR) = BBT12(theta_all(i), r_all(i)) / dr / d_theta;
        temp(UL) = -BBT12(theta_all(i), r_all(i)) / dr / d_theta;

        temp(UUR) = -BBT12(theta_all(i), r_all(i)) / dr / 4 / d_theta;
        temp(UUL) = BBT12(theta_all(i), r_all(i)) / dr / 4 / d_theta;


        if theta_all(i) == theta(1)
            b(i) = b(i) + T_bar * (BBT12(theta_all(i), r_all(i)) / dr / d_theta...
                                       -BBT12(theta_all(i), r_all(i)) / dr / 4 / d_theta);
        end

        if theta_all(i) == theta(end)
            b(i) = b(i) + T_bar * (BBT12(theta_all(i), r_all(i)) / dr / d_theta...
                                       -BBT12(theta_all(i), r_all(i)) / dr / 4 / d_theta);
        end

        temp_bc(i) = BBT22(theta_all(i), r_all(i)) / dr * (-3/2);
        temp_bc(right) = BBT12(theta_all(i), r_all(i)) / 2 / d_theta;
        temp_bc(left) = -BBT12(theta_all(i), r_all(i)) / 2 / d_theta;
        temp_bc(up) = BBT22(theta_all(i), r_all(i)) / dr * (2);
        temp_bc(up2) = BBT22(theta_all(i), r_all(i)) / dr * (-1/2);

        if theta_all(i) == theta(1)
            b(N+counter) = b(N+counter) + T_bar * (BBT12(theta_all(i), r_all(i)) / 2 / d_theta);
        end

        if theta_all(i) == theta(end)
            b(N+counter) = b(N+counter) + T_bar * (BBT12(theta_all(i), r_all(i)) / 2 / d_theta);
        end

        L_dagger(i, :) = temp;
        L_dagger(N+counter, :) = temp_bc;
        counter = counter + 1;
        disp(i/N*100)
    end
end

T = L_dagger \ b;
T = T - min(T);

figure 
plot3(theta_all, r_all, T, '.')
xlabel('\theta')
ylabel('r')
zlabel('T')
set(gca, 'FontSize', 20)
axis square

theta = reshape(theta_all, m, m);
r = reshape(r_all, m, m);
T = reshape(T, m, m);

X = r .* cos(theta);
Y = r .* sin(theta);

figure
contourf(X, Y, T, 60)
xlabel('x')
ylabel('y')
set(gca, 'FontSize', 20)
axis square
colorbar

% save('T_matrix', 'X', 'Y', 'T')
