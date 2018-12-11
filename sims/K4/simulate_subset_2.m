common;

rates = zeros(length(SNR), 1);

S = 2;

for i = 1:length(SNR)
    snr = SNR(i);

    rate = cacsim_subset(L, K, S, 2, snr, 'realizations', REALIZATIONS);

    [rel, iter] = size(rate);

    tmp = sum(rate, 1);
    rates(i) = rates(i) + (tmp(end) / rel);

    % Store the results 
    name = 'CC-BF-SCA (\alpha = \beta = 1)';
    save('data/sim-subset-2.mat', 'rates', 'name', 'SNR')
end
