common;

rates = zeros(length(SNR), 1);

for i = 1:length(SNR)
    snr = SNR(i)

    rate = cacsim_mu(L, K, snr, 'realizations', REALIZATIONS);

    [rel, iter] = size(rate);

    tmp = sum(rate, 1);
    rates(i) = rates(i) + (tmp(end) / rel);

    % Store the results 
    name = 'Unicast MaxMin SINR';
    save('data/sim-unicast.mat', 'rates', 'name', 'SNR')
end
