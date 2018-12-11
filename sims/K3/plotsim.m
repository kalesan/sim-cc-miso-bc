sims = char('data/sim-subset-3.mat', 'data/sim-subset-2-beta-1.mat', 'data/sim-subset-3-zf.mat', 'data/sim-subset-3-zf-eq.mat', ...
             'data/sim-unicast.mat');

common

figure,
hold on

altstyles = {'bo-'; 'bs--'; 'rd-'; 'r+--'; 'kv-'; 'k*--'; 'mv-'};

[n, ~] = size(sims);
for i = 1:n
    load(strtrim(sims(i, :)))

    plot(SNR, rates / log(2), altstyles{i}, 'LineWidth', 1.5, 'MarkerSize', 8);
end
grid on

legend('CC-SCA (\alpha = \beta = 2)', 'CC-SCA (\alpha = \beta = 1)', ...
        'CC-ZF (power loading)', ...
        'CC-ZF (equal power)', ...
         'Unicast MaxMin SINR', 'Location', 'northwest');
ylabel('Average rate per UE [bits/s/Hz]');
xlabel('SNR [dB]');
%title(['Interference coordination with caching ' ...  
       %'L = 2, K = ' int2str(K)]);

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
savefig(gcf, 'figs/K3.fig')
saveas(gcf, 'figs/K3.eps', 'epsc')
saveas(gcf, 'figs/K3.pdf')
saveas(gcf, 'figs/K3.png')
