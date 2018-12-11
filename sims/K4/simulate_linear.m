common;

rates = zeros(length(SNR), 1);

for i = 1:length(SNR)

    snr = SNR(i)

    rate = cacsim_subset_linear(L, K, snr, 'realizations', REALIZATIONS);

    [rel, iter] = size(rate);

    tmp = sum(rate, 1);
    rates(i) = rates(i) + (tmp(end) / rel);

    % Store the results 
    name = 'CC-BF-SCA (\alpha = 3, \beta = 1)';
    save('data/sim-alpha-3-beta-1.mat', 'rates', 'name', 'SNR')
end
