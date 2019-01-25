function rate = cacsim_subset_linear(L, K, SNR_dB, varargin)
    % Simulator for alpha = K-1 and beta=1, i.e., zero cross-user interference.

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

    MAI = cell(2, 1);
    MAI{1} = MA(1:(K-1), :);

    [N, ~] = size(MAI{1});

    tmp = zeros(N, K-2);

    for n = 1:N
        MA = MAI{1};
        tmp(n, :) = setdiff(1:K, MA(n, :));
    end

    MAI{2} = tmp;

    [msgs, ~] = size(MA);

    rate = zeros(REALIZATIONS, ITERATIONS);

    for (r = 1:REALIZATIONS)
        randn('seed', r);

        h = (1/sqrt(2)) * (randn(L, K) + 1j*randn(L, K));

        m_p = 1/sqrt(2) * (randn(L, msgs, 2) + randn(L, msgs, 2)*1j);

        m_p = m_p / (norm(m_p(:))) * sqrt(SNR);

        y_p = SINR(h, m_p, N0, MAI);

        R_prev = inf;
        
        for (i = 1:ITERATIONS)
            cvx_begin
                cvx_solver('sedumi');

                cvx_quiet('true')

                variable y(msgs, K) nonnegative;

                variable m(L, msgs, 2) complex;

                variable t(msgs, 1);

                maximize(sum(t(:)));

                subject to
                    for msg = 1:msgs
                        y(msg, :) >= t(msg);
                        sum_square_abs(vec(m(:, msg, :))) <= SNR;
                    end

                    % SINR Constraints
                    for msg = 1:msgs
                        MA = MAI{1};
                        for k = MA(msg, :)
                            c = N0 + norm(h(:, k)' * m_p(:, msg, 2))^2 + ...
                                norm(h(:, k)' * m_p(:, msg, 1))^2;

                            lhs = c / (1 + y_p(msg, k)) ...
                                  + (c / (1 + y_p(msg, k))^2) * (y_p(msg, k) - y(msg, k));

                            lhs = lhs - (2 / (1 + y_p(msg, k))) * ...
                                real(m_p(:, msg, 1)' * h(:, k) * h(:, k)' * (m_p(:, msg, 1) - m(:, msg, 1)));
                            lhs = lhs - (2 / (1 + y_p(msg, k))) * ...
                                real(m_p(:, msg, 2)' * h(:, k) * h(:, k)' * (m_p(:, msg, 2) - m(:, msg, 2)));

                            % Noise plus interference from group 2
                            rhs = N0;
                            rhs = rhs + sum_square_abs(h(:, k)' * m(:, msg, 2));

                            lhs >= rhs;
                        end

                        MA = MAI{2};
                        for k = MA(msg, :)
                            c = N0 + norm(h(:, k)' * m_p(:, msg, 2))^2 + ...
                                norm(h(:, k)' * m_p(:, msg, 1))^2;

                            lhs = c / (1 + y_p(msg, k)) ...
                                  + (c / (1 + y_p(msg, k))^2) * (y_p(msg, k) - y(msg, k));

                            lhs = lhs - (2 / (1 + y_p(msg, k))) * ...
                                real(m_p(:, msg, 2)' * h(:, k) * h(:, k)' * (m_p(:, msg, 2) - m(:, msg, 2)));

                            lhs = lhs - (2 / (1 + y_p(msg, k))) * ...
                                real(m_p(:, msg, 1)' * h(:, k) * h(:, k)' * (m_p(:, msg, 1) - m(:, msg, 1)));

                            % Noise plus interference from group 1
                            rhs = N0;
                            rhs = rhs + sum_square_abs(h(:, k)' * m(:, msg, 1));

                            lhs >= rhs;
                        end
                    end, 
            cvx_end

            m_p = m;

            % SINR
            y_p = SINR(h, m_p, N0, MAI);

            M = sum(MA(:) == 1);

            % Rate for each message
            T_m = (1/K) ./ (log(1 + min(y_p'))); 

            R = 1 / sum(T_m);

            rate(r, i) = R;

            disp(['Realization ' int2str(r) ' rate ' num2str(rate(r, i))  ...
                  ' pwr ' num2str(norm(m_p(:))^2 / (msgs * SNR) * 100) '%'])

            if (abs(rate(r, i) - R_prev) < 1e-4)
                rate(r, i:end) = rate(r, i);

                break;
            end

            R_prev = rate(r, i);
        end
    end
end

function y = SINR(h, m, N0, MAI)
    [L, K] = size(h);
    [~, msgs, T] = size(m);

    y = inf*ones(msgs, K);

    for msg = 1:msgs;
        MA = MAI{1};

        for k = MA(msg, :)
            y(msg, k) = norm(h(:, k)' * m(:, msg, 1))^2 / ...
                        (N0 + norm(h(:, k)' * m(:, msg, 2))^2);
        end

        MA = MAI{2};
        for k = MA(msg, :)
            y(msg, k) = norm(h(:, k)' * m(:, msg, 2))^2 / ...
                        (N0 + norm(h(:, k)' * m(:, msg, 1))^2);
        end
    end
end
