common;

rates = zeros(length(SNR), 1);

% Subset size (alpha + 1)
SS = 3;

for i = 1:length(SNR)
    snr = SNR(i);

    rate = cacsim_subset(L, K, SS, S, snr, 'realizations', realizations);

    [rel, iter] = size(rate);

    tmp = sum(rate, 1);
    rates(i) = rates(i) + (tmp(end) / rel);

    % Store the results 
    name = '\alpha = \beta = 2';
    save('data/sim-subset-3.mat', 'rates', 'name', 'SNR')
end
