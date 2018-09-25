clear

T_bar = 16.233624025276129; % from flux

figure

coordinates = load('isophase_init_conds.txt');
r = sqrt(coordinates(:,1).^2 + coordinates(:,2).^2);

plot([min(r)-0.1 max(r)+0.1], [T_bar T_bar], 'r-', 'LineWidth', 5)
hold on

load isophase_data.mat

plot(r, means_all, 'k.-', 'MarkerSize', 45)
% errorbar(r, means_all, se_all,  'k', 'LineWidth', 10)

coordinates = load('spoke_init_conds.txt');
r = sqrt(coordinates(:,1).^2 + coordinates(:,2).^2);
load spoke_data.mat

plot(r, means_all, 'ksq-', 'MarkerSize', 15, 'MarkerFaceColor', 'k')
% errorbar(r, means_all, se_all,  'k', 'LineWidth', 10)
axis([min(r)-0.1 max(r)+0.1 8 17])
axis square
grid on
x = xlabel('Radial distance from origin r');
y = ylabel('Mean first-return time');
set(x,'Interpreter','latex','fontsize',20)
set(y,'Interpreter','latex','fontsize',20)
set(gca, 'FontSize', 20)

h = legend('Mean period $$\overline{T}$$', 'Isophase', 'Spoke', 'location', 'SouthEast');
set(h,'Interpreter','latex','fontsize',15)

print('-depsc', 'het-compare')

