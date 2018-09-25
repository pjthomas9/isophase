clear

data = load('clock_model_traj.dat');
x = data(:,4);
y = data(:,5);

figure
plot(x,y, 'k')
x = xlabel('x');
y = ylabel('y');
set(x,'Interpreter','latex','fontsize',20)
set(y,'Interpreter','latex','fontsize',20)
axis([-1.5 1.5 -1.5 1.5])
set(gca, 'FontSize', 20)
axis square
grid on
print('-depsc', 'clock-model-traj')