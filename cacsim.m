function rate = cacsim(L, K, SNR_dB, varargin)
    p = inputParser;
    
    p.addOptional('realizations',  500, @isnumeric);
    p.addOptional('iterations',  100, @isnumeric);

    p.parse(varargin{:});
    params = p.Results;

    REALIZATIONS = params.realizations;
    ITERATIONS   = params.iterations;

    SNR = 10^(SNR_dB / 10);

    N0 = 1;

    if (K ~= 3)
        error('Current implementation supports only 3 UE');
    end

    rate = zeros(REALIZATIONS, ITERATIONS);

    for (r = 1:REALIZATIONS)
        randn('seed', r);

        h = (1/sqrt(2)) * (randn(L, K) + randn(L, K));

        m_p = 1/sqrt(2) * (randn(L, K) + randn(L, K)*1j);
        m_p = m_p / norm(m_p(:)) * sqrt(SNR);

        y_p = SINR(h, m_p, N0);

        R_prev = inf;

        for (i = 1:ITERATIONS)
            cvx_begin
                cvx_solver('sedumi');

                cvx_quiet('true')

                variable y(K, 2) nonnegative;

                variable m(L, K) complex;

                variable eph;

                maximize(eph);

                subject to
                    % Rate bounds
                    log(1 + sum(y, 2)) / log(2) >= eph;
                    log(1 + y(:, 1)) / log(2) >= eph;
                    log(1 + y(:, 2)) / log(2) >= eph;

                    % Power constraint
                    sum_square_abs(m(:)) <= SNR;

                    % SINR Constraints
                    for k = 1:K; t = 1; for j = [1:(k-1), (k+1):K]
                        c = N0 + norm(h(:, k)' * m_p(:, j))^2 + ...
                                    norm(h(:, k)' * m_p(:, k))^2;

                        c  - 2 * real(m_p(:, j)' * h(:, k) * h(:, k)' * (m_p(:, j) - m(:, j))) ...
                            - 2 * real(m_p(:, k)' * h(:, k) * h(:, k)' * (m_p(:, k) - m(:, k))) ...
                            + (c / (1 + y_p(k, t))) * (y_p(k, t) - y(k, t)) ...
                                >= (1 + y_p(k ,t)) * (sum_square_abs(h(:, k)' * m(:, k)) + N0);

                        t = t + 1;
                    end, end
            cvx_end

            m_p = m;

            % SINR
            y_p = SINR(h, m, N0);

            R_sum = 0.5 * log2(1 + sum(y_p, 2));
            R_1 = log2(1 + y_p(:, 1));
            R_2 = log2(1 + y_p(:, 2));

            R = min([R_sum(:); R_1(:); R_2(:)]);

            rate(r, i) = R;

            disp(['Realization ' int2str(r) ' rate ' num2str(rate(r, i))])

            if (abs(rate(r, i) - R_prev) < 1e-4)
                rate(r, i:end) = rate(r, i);

                break;
            end

            R_prev = rate(r, i);
        end
    end
end

function y = SINR(h, m, N0)
    [L, K] = size(h);

    y = zeros(K, 2);
    for k = 1:K; 
        t = 1; 
        
        for j = [1:(k-1), (k+1):K]
            y(k, t) = norm(h(:, k)' * m(:, j))^2 / ...
                        (N0 + norm(h(:, k)' * m(:, k))^2);

            t = t + 1;
        end
    end
end
