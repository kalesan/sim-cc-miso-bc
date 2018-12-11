function rate = cacsim_subset_K6B3(L, K, SS, S, SNR_dB, varargin)
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

    groups = K / beta;

    K_g = K / groups;

    MAI = nchoosek(1:K_g, S);

    groupings = nchoosek(K, K_g) / 2;
    
    % Take generalized caching for subset of size 
    MAI = nchoosek(1:K_g, S);

    [msgs, ~] = size(MAI);
    [uem] = sum(MAI(:, 1) == 1);

    alloc = nchoosek(1:K, SS)';

    rate = zeros(REALIZATIONS, ITERATIONS);

    % Transmission time slots
    T = nchoosek(K, SS);

    for (r = 1:REALIZATIONS)
        randn('seed', r);

        h = (1/sqrt(2)) * (randn(L, K) + 1j*randn(L, K));

        m_p = 1/sqrt(2) * (randn(L, msgs, groups, groupings) + randn(L, msgs, groups, groupings)*1j);

        for t = 1:groupings
            tmp = m_p(:, :, :, t);
            m_p(:, :, :, t) = m_p(:, :, :, t) / (norm(tmp(:))) * sqrt(SNR);
        end

        y_p = SINR(h, m_p, N0, alloc, SS, S);

        R_prev = inf;

        for (i = 1:ITERATIONS)
            m = zeros(size(m_p));

            m = solve_socp_beta_K6B3(h, m_p, y_p, SNR, N0, S, beta);

            m_p = m;

            % SINR
            y_p = SINR(h, m, N0, alloc, SS, S);

            T = zeros(groupings, 1);
            for g = 1:groupings
                % MAC rate points R_MAC
                R_mac = [];
                for k = 1:K_g, for v = nchoosek(1:uem, k)'
                    tmp = log(1 + sum(y_p(:, v, :, g), 2)) / k;
                    R_mac = [R_mac tmp];
                end, end

                T(g) = T(g) + (1 / min(R_mac(:)));
            end

            alpha = 5;
            beta = 2;
            t = 1;

            gamma = (alpha+1) / (beta + 1);

            n = nchoosek(K, t) * nchoosek(K-t-1, alpha-1) * factorial(alpha-1) / ...
                (factorial(gamma-1) * factorial(beta-1) * factorial(beta+1)^(gamma-1));

            rate(r, i) = n/sum(T(:));

            disp(['Realization ' int2str(r) ' rate ' num2str(rate(r, i))  ...
                  ' pwr ' num2str(norm(m_p(:))^2 / (groupings * SNR) * 100) '%'])

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
    [~, msgs, groups, groupings] = size(m);

    K_g = K / groups;

    G = zeros(K_g, groupings, groups);
    tmp = nchoosek(1:K, K_g)';
    G(:, :, 1) = tmp(:, 1:groupings); 

    for g = 1:groupings
        G(:, g, 2) = setdiff(1:K, G(:,g,1));
    end

    MAI = nchoosek(1:K_g, S);

    [uem] = sum(MAI(:, 1) == 1);
    y = zeros(K_g, uem, groups, groupings);

    for g = 1:groupings, for f = 1:groups
        for k = 1:K_g
            t = 1;

            h_ = h(:, G(k, g, f));

            for msg = find(any(MAI' == k))
                % Interfering messages indices for k
                IM = find(all(MAI' ~= k));

                a = norm(h_' * m(:, msg, f, g))^2;
                b = N0 + sum_square_abs(h_' * m(:, IM, f, g));;

                for f_ = [1:(f-1),(f+1):groups]
                    b = b + sum_square_abs(vec(h_' * m(:, :, f_, g)));
                end

                y(k, t, f, g) = a / b;

                t = t + 1;
            end
        end
    end, end
end
