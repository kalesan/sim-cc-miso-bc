function rate = cc(precoder, filename, L, K, S, SNR_dB, varargin)
    % Generic coded caching simulator that uses precoder type objects for
    % specific precoder design.

    p = inputParser;
    
    p.addOptional('realizations',  500, @isnumeric);
    p.addOptional('iterations',  100, @isnumeric);
    p.addOptional('name', '', @ischar);

    p.parse(varargin{:});
    params = p.Results;

    REALIZATIONS = params.realizations;
    ITERATIONS   = params.iterations;

    if length(params.name) == 0
        name = filename;
    else
        name = params.name;
    end

    N0 = 1;

    MAI = nchoosek(1:K, S);

    [msgs, ~] = size(MAI);
    [uem] = sum(MAI(:, 1) == 1);

    rates = zeros(length(SNR_dB), 1);

    for s = 1:length(SNR_dB)
        SNR = 10^(SNR_dB(s) / 10);

        rate = zeros(REALIZATIONS, ITERATIONS);

        for (r = 1:REALIZATIONS)
            randn('seed', r);

            h = (1/sqrt(2)) * (randn(L, K) + 1j*randn(L, K));

            disp(['Realization ' int2str(r) '/' int2str(REALIZATIONS)])

            rate(r, :) = precoder.simulate(h, S, SNR, N0, MAI, ITERATIONS);
        end

        [rel, iter] = size(rate);

        tmp = sum(rate, 1);
        rates(s) = rates(s) + (tmp(end) / rel);
        
        % Store the results 
        save(filename, 'rates', 'name', 'SNR_dB')
    end
end
