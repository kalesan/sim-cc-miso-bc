common;

rates = zeros(length(SNR), 1);

% Subset size (alpha + 1 = 6)
SS = K;

for i = 1:length(SNR)
    snr = SNR(i);

    rate = cacsim_subset_K6B2_zf(L, K, SS, S, snr, 'realizations', realizations, 'beta', 2);

    [rel, iter] = size(rate);

    tmp = sum(rate, 1);
    rates(i) = rates(i) + (tmp(end) / rel);

    % Store the results 
    name = 'CC-ZF (\alpha = 5, \beta = 1)';
    save('data/sim-subset-6-beta-1-zf.mat', 'rates', 'name', 'SNR')
end
