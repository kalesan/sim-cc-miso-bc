sims = char('data/sim-subset-6.mat', 'data/sim-subset-5.mat', 'data/sim-subset-4.mat', ...
            'data/sim-subset-3.mat', 'data/sim-subset-2.mat');

common

SNR = -10:5:40;

figure,

altstyles = {'bo-'; 'rd-'; 'k+-'; 'mv-'; 'g*-'; 'mv--'; 'bs--'};

[n, ~] = size(sims);
hold on
for i = 1:n
    load(strtrim(sims(i, :)))
    SNR = -10:5:40;

    plot(SNR, rates(1:length(SNR)) / log(2), altstyles{i}, 'LineWidth', 1.5, 'MarkerSize', 8);
end
hold off
grid on
legend('\alpha = \beta = 5', '\alpha = \beta = 4', '\alpha = \beta = 3',  ...
        '\alpha = \beta = 2','\alpha = \beta = 1', 'Location', 'northwest') 
ylabel('Average rate per UE [bits/s/Hz]');
xlabel('SNR [dB]');

% Font and line sizes
set(gca, 'FontSize', 12)
set(gca, 'LineWidth', 1.5)

% Fill the figure
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% Save the figure
savefig(gcf, 'figs/K6-alpha.fig')
saveas(gcf, 'figs/K6-alpha.eps', 'epsc')
saveas(gcf, 'figs/K6-alpha.pdf')
saveas(gcf, 'figs/K6-alpha.png')
