sims = char('data/sim-subset-3.mat', 'data/sim-subset-3-l2.mat', 'data/sim-subset-3-zf.mat', ...
            'data/sim-unicast.mat');

common

figure,
hold on

altstyles = {'bo-'; 'bs--'; 'rd-'; 'kv-'; 'r+--'; 'k*--'; 'mv-'};

[n, ~] = size(sims);
for i = 1:n
    load(strtrim(sims(i, :)))

    plot(SNR, rates / log(2), altstyles{i}, 'LineWidth', 1.5, 'MarkerSize', 8);
end
grid on

legend('CC-SCA (\alpha = \beta = 2)', 'CC-SCA (\alpha = \beta = 2) [L=2]', ...
       'CC-ZF (power loading)', ...
       'Unicast MaxMin SINR', 'Location', 'northwest');
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
savefig(gcf, 'figs/K3L3.fig')
saveas(gcf, 'figs/K3L3.eps', 'epsc')
saveas(gcf, 'figs/K3L3.pdf')
saveas(gcf, 'figs/K3L3.png')
