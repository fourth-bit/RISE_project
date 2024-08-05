% Graphing Waveforms

figure
space = linspace(0, 10, 10001);
plot(space, cos(space * 2 * pi), 'LineWidth', 3)
ylim([-1 2])