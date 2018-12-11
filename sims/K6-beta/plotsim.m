sims = char('data/sim-subset-6.mat', 'data/sim-subset-6-beta-1.mat', 'data/sim-subset-6-beta-2.mat', ...
            'data/sim-subset-6-zf.mat', 'data/sim-subset-6-beta-1-zf.mat', 'data/sim-subset-6-beta-2-zf.mat');

common

SNR = -10:5:40;

figure,

altstyles = {'bo-'; 'rd-'; 'k+-'; 'bo--'; 'rd--'; 'k+--'; 'bs--'};

[n, ~] = size(sims);
hold on
for i = 1:n
    load(strtrim(sims(i, :)))
    SNR = -10:5:40;

    plot(SNR, rates(1:length(SNR)) / log(2), altstyles{i}, 'LineWidth', 1.5, 'MarkerSize', 8);
end
hold off
grid on
legend('CC-SCA (\alpha = \beta = 5)', 'CC-SCA (\alpha = 5, \beta = 1)', 'CC-SCA (\alpha = 5, \beta = 2)', ...
    'CC-ZF (\alpha = \beta = 5)', 'CC-ZF (\alpha = 5, \beta = 1)', 'CC-ZF (\alpha = 5, \beta = 2)', ...
        'Location', 'northwest');
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
savefig(gcf, 'figs/K6-beta.fig')
saveas(gcf, 'figs/K6-beta.eps', 'epsc')
saveas(gcf, 'figs/K6-beta.pdf')
saveas(gcf, 'figs/K6-beta.png')
