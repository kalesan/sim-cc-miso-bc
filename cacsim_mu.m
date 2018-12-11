function rate = cacsim_mu(L, K, SNR_dB, varargin)
    p = inputParser;
    
    p.addOptional('realizations',  500, @isnumeric);
    p.addOptional('iterations',  100, @isnumeric);

    p.parse(varargin{:});
    params = p.Results;

    REALIZATIONS = params.realizations;
    ITERATIONS   = params.iterations;

    SNR = 10^(SNR_dB / 10);

    N0 = 1;

    [MA, IA] = cache(K);

    [msgs, ~] = size(MA);
    [uem] = sum(MA(:, 1) == 1);

    rate = zeros(REALIZATIONS, ITERATIONS);

    for (r = 1:REALIZATIONS)
        randn('seed', r);

        h = (1/sqrt(2)) * (randn(L, K) + 1j*randn(L, K));

        m_p = 1/sqrt(2) * (randn(L, K) + randn(L, K)*1j);
        m_p = m_p / norm(m_p(:)) * sqrt(SNR); 
        y_p = SINR(h, m_p, N0);

        R_prev = inf;

        for (iter = 1:ITERATIONS)
            cvx_begin
                cvx_solver('sedumi');

                cvx_quiet('true')

                variable y(K, 1) nonnegative;

                variable m(L, K) complex;

                variable eph;

                maximize(eph);

                subject to
                    % Rate bounds
                    y >= eph;

                    % Power constraint
                    sum_square_abs(m(:)) <= SNR;

                    % v is the set of active UEs 
                    for k = 1:K
                        c = N0;
                        % j goes through all active users
                        for j = 1:K
                            c = c + norm(h(:, k)' * m_p(:, j))^2;
                        end

                        lhs = c / (1 + y_p(k)) ...
                              + (c / (1 + y_p(k))^2) * (y_p(k) - y(k));

                        % j goes through all active users
                        for j = 1:K
                            lhs = lhs - (2 / (1 + y_p(k))) * ...
                                real(m_p(:, j)' * h(:, k) * h(:, k)' * (m_p(:, j) - m(:, j)));
                        end

                        rhs = N0;

                        % i are the interfering users
                        for i = 1:K
                            if (k == i)
                                continue;
                            end

                            rhs = rhs + sum_square_abs(h(:, k)' * m(:, i));
                        end

                        lhs >= rhs;
                    end
            cvx_end

            m_p = m;

            % SINR
            y_p = SINR(h, m, N0);

            % Messages per UE
            M = sum(MA(:) == 1);

            R = ((M + 1) / M) * min(log(1 + y_p(:)));

            rate(r, iter) = R;

            disp(['Realization ' int2str(r) ' rate ' num2str(rate(r, iter))])

            if (abs(rate(r, iter) - R_prev) < 1e-4)
                rate(r, iter:end) = rate(r, iter);

                break;
            end

            R_prev = rate(r, iter);
        end
    end
end

function y = SINR(h, m, N0)
    [L, K] = size(h);

    y = inf*ones(K, 1);

    for k = 1:K
        I = N0;

        % j goes through the interfering UEs
        for j = 1:K
            if (j ~= k)
                I = I + norm(h(:, k)' * m(:, j))^2;
            end
        end

        y(k) = norm(h(:, k)' * m(:, k))^2 / I;
    end
end
