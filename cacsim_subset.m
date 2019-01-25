function rate = cacsim_subset(L, K, SS, S, SNR_dB, varargin)
    % Simulator for alpha=beta=SS scenarios. 
    %
    % Note that the simultaneously served user subset size is generally alpha+1.

    p = inputParser;
    
    p.addOptional('realizations',  500, @isnumeric);
    p.addOptional('iterations',  100, @isnumeric);
    p.addOptional('beta', [], @isnumeric);

    p.parse(varargin{:});
    params = p.Results;

    REALIZATIONS = params.realizations;
    ITERATIONS   = params.iterations;
    beta         = params.beta;

    if isempty(beta)
        beta = SS;
    end

    SNR = 10^(SNR_dB / 10);

    N0 = 1;

    % Take generalized caching for subset of size 
    MAI = nchoosek(1:SS, S);

    [msgs, ~] = size(MAI);
    [uem] = sum(MAI(:, 1) == 1);

    alloc = nchoosek(1:K, SS)';

    rate = zeros(REALIZATIONS, ITERATIONS);

    % Transmission time slots
    T = nchoosek(K, SS);

    for (r = 1:REALIZATIONS)
        randn('seed', r);

        h = (1/sqrt(2)) * (randn(L, K) + 1j*randn(L, K));

        m_p = 1/sqrt(2) * (randn(L, msgs, T) + randn(L, msgs, T)*1j);

        for t = 1:T
            tmp = m_p(:, :, t);
            m_p(:, :, t) = m_p(:, :, t) / (norm(tmp(:))) * sqrt(SNR);
        end

        y_p = SINR(h, m_p, N0, alloc, SS, S);

        R_prev = inf;

        for (i = 1:ITERATIONS)
            m = zeros(size(m_p));

            for t = 1:T
                h_ = zeros(L, SS);
                for k = 1:SS
                    h_(:, k) = h(:, alloc(k, t));
                end

                if (K == 6 && beta == 3)
                    m(:, :, t) = solve_socp_beta_K6B3(h_, m_p(:, :, t), y_p(:, :, t), SNR, N0, S, beta);
                else
                    m(:, :, t) = solve_socp(h_, m_p(:, :, t), y_p(:, :, t), SNR, N0, S);
                end
            end

            m_p = m;

            % SINR
            y_p = SINR(h, m, N0, alloc, SS, S);

            R = 0;
            for t = 1:T
                % MAC rate points R_MAC
                R_mac = [];
                if uem == 1
                    tmp = log(1 + y_p(:, :, t));
                    R_mac = [R_mac tmp(:)];
                else
                    for k = 1:SS, for v = nchoosek(1:uem, k)'
                        tmp = log(1 + sum(y_p(:, v, t), 2)) / k;
                        R_mac = [R_mac tmp];
                    end, end
                end

                R = R + (1 / min(R_mac(:)));
            end

            split = nchoosek(K, SS) * max(1, nchoosek(K - 2, SS - 2));

            % Number of times we are transmitting with distinct user pair in
            % place.
            subfiles = nchoosek(K - 2, SS - 2);

            split = nchoosek(K, SS) * max(1, split);

            % Total number of messages (per file). Each user has subfiles
            % number of parts.
            split = subfiles * K;

            rate(r, i) = split/R;

            disp(['Realization ' int2str(r) ' rate ' num2str(rate(r, i))  ...
                  ' pwr ' num2str(norm(m_p(:))^2 / (T * SNR) * 100) '%'])

            if (abs(rate(r, i) - R_prev) < 1e-4)
                rate(r, i:end) = rate(r, i);

                break;
            end

            R_prev = rate(r, i);
        end
    end
end

function y = SINR(h, m, N0, alloc, SS, S)
    [L, K] = size(h);
    [~, msgs, T] = size(m);

    MAI = nchoosek(1:SS, S);

    [uem] = sum(MAI(:, 1) == 1);
    y = zeros(SS, uem, T);

    for slot = 1:T
        h_ = zeros(L, SS);
        for k = 1:SS
            h_(:, k) = h(:, alloc(k, slot));
        end

        y(:, :, slot) = sinr(h_, m(:, :, slot), N0, S);
    end
end
