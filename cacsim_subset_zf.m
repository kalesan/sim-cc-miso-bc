function rate = cacsim_subset_zf(L, K, S, SNR_dB, varargin)
    p = inputParser;
    
    p.addOptional('realizations',  500, @isnumeric);
    p.addOptional('iterations',  100, @isnumeric);

    p.parse(varargin{:});
    params = p.Results;

    REALIZATIONS = params.realizations;
    ITERATIONS   = params.iterations;

    SNR = 10^(SNR_dB / 10);

    N0 = 1;

    % Take generalized caching for subset of size 
    [MA, IA] = cache(S);

    [msgs, ~] = size(MA);
    [uem] = sum(MA(:, 1) == 1);

    % In ZF we need to map the msgs to actual users
    [MA_zf, IA_zf] = cache(K);
    [msgs_zf, ~] = size(MA_zf);
    [uem_zf] = sum(MA_zf(:, 1) == 1);

    alloc = nchoosek(1:K, S)';

    rate = zeros(REALIZATIONS, ITERATIONS);

    % Transmission time slots
    T = nchoosek(K, S);

    msg_index = zeros(nchoosek(S, 2), T);
    for t = 1:T
        tmp = nchoosek(alloc(:, t), 2);

        for k = 1:nchoosek(S, 2)
            msg_index(k, t) = find(ismember(MA_zf, tmp(k, :), 'rows'));
        end
    end

    for (r = 1:REALIZATIONS)
        randn('seed', r);

        h = (1/sqrt(2)) * (randn(L, K) + 1j*randn(L, K));

        m_zf = zeros(L, K);

        for k = 1:msgs_zf
            tmp = null(h(:, any(IA_zf == k, 2))');
            [U, D, V] = svd(h'  * tmp);

            m_zf(:, k) = tmp  * V(:, 1);
            m_zf(:, k) = m_zf(:, k) / norm(m_zf(:, k));
        end

        m_p = 1/sqrt(2) * (randn(L, msgs, T) + randn(L, msgs, T)*1j);

        for t = 1:T
            tmp = m_p(:, :, t);
            m_p(:, :, t) = m_zf(:, msg_index(:, t)) ...
                            / (norm(tmp(:))) * sqrt(SNR / msgs);
        end

        y_p = SINR(h, m_p, N0, MA, IA, alloc, S);

        R_prev = inf;

        for (i = 1:ITERATIONS)
            m = zeros(size(m_p));

            for t = 1:T
                h_ = zeros(L, S);
                for k = 1:S
                    h_(:, k) = h(:, alloc(k, t));
                end

                m(:, :, t) = solve_zf(h_, MA, IA, m_p(:, :, t), ...
                                      m_zf(:, msg_index(:, t)), ...
                                      y_p(:, :, t), SNR, N0);
            end

            m_p = m;

            % SINR
            y_p = SINR(h, m, N0, MA, IA, alloc, S);

            R = 0;
            for t = 1:T
                % MAC rate points R_MAC
                R_mac = [];
                for k = 1:S, for v = nchoosek(1:uem, k)'
                    tmp = log(1 + sum(y_p(:, v, t), 2)) / k;
                    R_mac = [R_mac tmp];
                end, end

                R = R + (1 / min(R_mac(:)));
            end

            split = nchoosek(K, S) * max(1, nchoosek(K - 2, S - 2));

            % Number of times we are transmitting with distinct user pair in
            % place.
            subfiles = nchoosek(K - 2, S - 2);

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

function y = SINR(h, m, N0, MA, IA, alloc, S)
    [L, K] = size(h);
    [~, msgs, T] = size(m);

    [uem] = sum(MA(:, 1) == 1);
    y = zeros(S, uem, T);

    for slot = 1:T
        h_ = zeros(L, S);
        for k = 1:S
            h_(:, k) = h(:, alloc(k, slot));
        end

        y(:, :, slot) = sinr(h_, m(:, :, slot), N0, MA, IA);
    end
end
