common;

rates = zeros(length(SNR), 1);

for i = 1:length(SNR)

    snr = SNR(i)

    rate = cacsim_zf_equal_power(L, K, snr, 'realizations', REALIZATIONS);

    rates(i) = rates(i) + mean(rate(:,end));

    % Store the results 
    name = 'CC-ZF (equal power)';
    save('data/sim-subset-3-zf-eq.mat', 'rates', 'name', 'SNR')
end

