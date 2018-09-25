clear

T_bar = 20.929154018932547; % from flux

figure

coordinates = load('isophase_init_conds.txt');
r = coordinates(:,2);

plot([min(r)-0.1 max(r)+0.1], [T_bar T_bar], 'r-', 'LineWidth', 5)
hold on

load isophase_data.mat

plot(r, means_all, 'k.-', 'MarkerSize', 45)
% errorbar(r, means_all, se_all,  'k', 'LineWidth', 10)

coordinates = load('spoke_init_conds.txt');
r = coordinates(:,2);
load spoke_data.mat

plot(r, means_all, 'ksq-', 'MarkerSize', 15, 'MarkerFaceColor', 'k')
% errorbar(r, means_all, se_all,  'k', 'LineWidth', 10)
axis([min(r)-0.1 max(r)+0.1 18 22])
axis square
grid on
x = xlabel('Radial distance from origin r');
y = ylabel('Mean first-return time');
set(x,'Interpreter','latex','fontsize',20)
set(y,'Interpreter','latex','fontsize',20)
set(gca, 'FontSize', 20)

h = legend('Mean period $$\overline{T}$$', 'Isophase', 'Spoke', 'location', 'SouthEast');
set(h,'Interpreter','latex','fontsize',15)

% print('-depsc', 'y-polar-SL-NE-compare')

