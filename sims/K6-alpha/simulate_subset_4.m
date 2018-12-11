common;

rates = zeros(length(SNR), 1);

% Subset size (alpha + 1)
SS = 4;

for i = 1:length(SNR)
    snr = SNR(i);

    rate = cacsim_subset(L, K, SS, S, snr, 'realizations', realizations);

    [rel, iter] = size(rate);

    tmp = sum(rate, 1);
    rates(i) = rates(i) + (tmp(end) / rel);

    % Store the results 
    name = '\alpha = \beta = 3';
    save('data/sim-subset-4.mat', 'rates', 'name', 'SNR')
end
