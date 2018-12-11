function rate = cacsim_subset_K6B2_zf(L, K, SS, S, SNR_dB, varargin)
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

    G = [];

    G1 = nchoosek(1:K, 2);
    for i = 1:(nchoosek(K, 2) / 3)
        a1 = G1(i, :);

        G2 = nchoosek(setdiff(1:K, a1), 2);
        for j = 1:(nchoosek(K-2, 2) / 2)
            a2 = G2(j, :);
            a3 = setdiff(1:K, [a1 a2]);

            G = [G; a1 a2 a3];
        end
    end

    G_ = G;

    [groupings, ~] = size(G);
    G = reshape(G_', [2, 3, groupings]);
    G = permute(G, [1 3 2]); % [K_g, groupings, groups]
    
    % Take generalized caching for subset of size 
    MAI = nchoosek(1:K_g, S);

    [msgs, ~] = size(MAI);

    alloc = nchoosek(1:K, SS)';

    rate = zeros(REALIZATIONS, ITERATIONS);

    % Transmission time slots
    T = nchoosek(K, SS);

    for (r = 1:REALIZATIONS)
        randn('seed', r);

        h = (1/sqrt(2)) * (randn(L, K) + 1j*randn(L, K));

        m_p = 1/sqrt(2) * (randn(L, groups, groupings) + randn(L, groups, groupings)*1j);

        for t = 1:groupings
            m_zf = zeros(L, groups);

            for f = 1:groups
                tmp = null(h(:, setdiff(1:K, G(:, t, f))')');
                [U, D, V] = svd(h'  * tmp);

                m_zf(:, f) = tmp  * V(:, 1);
                m_zf(:, f) = (m_zf(:, f) / norm(m_zf(:, f))) * (sqrt(SNR / groups));
            end

            m_p(:, :, t) = m_zf;
        end

        y_p = SINR(h, m_p, N0, alloc, SS, S, G);

        R_prev = inf;

        for (i = 1:ITERATIONS)
            m = zeros(size(m_p));

            m = solve_socp_beta_K6B2_zf(h, m_p, y_p, SNR, N0, S, beta, G);

            m_p = m;

            % SINR
            y_p = SINR(h, m, N0, alloc, SS, S, G);

            T = zeros(groupings, 1);
            for t = 1:groupings
                % MAC rate points R_MAC
                R_mac = [];

                tmp = log(1 + vec(y_p(:, :, t)));
                R_mac = [R_mac tmp(:)];

                T(t) = T(t) + (1 / min(R_mac(:)));
            end

            alpha = 5;
            beta = 1;
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

function y = SINR(h, m, N0, alloc, SS, S, G)
    [L, K] = size(h);
    [~, groups, groupings] = size(m);

    K_g = K / groups;

    MAI = nchoosek(1:K_g, S);

    y = zeros(K_g, groups, groupings);

    for g = 1:groupings, for f = 1:groups
        for k = 1:K_g
            h_ = h(:, G(k, g, f));

            % Interfering messages indices for k
            a = norm(h_' * m(:, f, g))^2;

            b = N0 + sum_square_abs(vec(h_' * m(:, [1:(f-1),(f+1):groups], g)));

            y(k, f, g) = a / b;
        end
    end, end
end
