common;

rates = zeros(length(SNR), 1);

for i = 1:length(SNR)

    snr = SNR(i)

    rate = cacsim_subset(L, K, K, S, snr, 'realizations', realizations);

    [rel, iter] = size(rate);

    tmp = sum(rate, 1);
    rates(i) = rates(i) + (tmp(end) / rel);

    % Store the results 
    name = 'CC-BF-SCA (\alpha = \beta = 2)';
    save('data/sim-subset-3.mat', 'rates', 'name', 'SNR')
end
