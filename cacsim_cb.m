function rate = cacsim_cb(L, K, SNR_dB, varargin)
    p = inputParser;
    
    p.addOptional('realizations',  500, @isnumeric);
    p.addOptional('iterations',  100, @isnumeric);

    p.parse(varargin{:});
    params = p.Results;

    REALIZATIONS = params.realizations;
    ITERATIONS   = params.iterations;

    SNR = 10^(SNR_dB / 10);

    N0 = 1;

    % Number of transmission time slots for each UE
    T = nchoosek(K, L);

    rate = zeros(REALIZATIONS, ITERATIONS);

    for (r = 1:REALIZATIONS)
        randn('seed', r);

        h = (1/sqrt(2)) * (randn(L, K) + 1j*randn(L, K));

        m_p = 1/sqrt(2) * (randn(L, K, T) + randn(L, K, T)*1j);
        m_p = m_p / norm(m_p(:)) * sqrt(SNR); 
        y_p = SINR(h, m_p, N0);

        R_prev = inf;

        for (iter = 1:ITERATIONS)
            cvx_begin
                cvx_solver('sedumi');

                cvx_quiet('true')

                variable y(K, T) nonnegative;

                variable m(L, K, T) complex;

                variable eph(T, 1);

                maximize(sum(eph));

                subject to
                    % Rate bounds
                    ts = 1;
                    for v = nchoosek(1:K, L)'
                        y(:, ts) >= eph(ts);
                        ts = ts + 1;
                    end

                    % Power constraint
                    for t = 1:T
                        sum_square_abs(vec(m(:, :, t))) <= SNR;
                    end

                    % Active time slot index
                    ts = 1;

                    % v is the set of active UEs 
                    for v = nchoosek(1:K, L)'
                        % Null out the non-active UEs 
                        for j = setdiff(1:K, v)
                            m(:, j, ts) == 0;
                        end

                        % k is the intended user 
                        for k = v'
                            c = N0;
                            % j goes through all active users
                            for j = v'
                                c = c + norm(h(:, k)' * m_p(:, j, ts))^2;
                            end

                            lhs = c / (1 + y_p(k , ts)) ...
                                  + (c / (1 + y_p(k, ts))^2) * (y_p(k, ts) - y(k, ts));

                            % j goes through all active users
                            for j = v'
                                lhs = lhs - (2 / (1 + y_p(k, ts))) * ...
                                    real(m_p(:, j, ts)' * h(:, k) * h(:, k)' * (m_p(:, j, ts) - m(:, j, ts)));
                            end
                            rhs = N0;

                            % i are the interfering users
                            for i = v'
                                if (k == i)
                                    continue;
                                end

                                rhs = rhs + sum_square_abs(h(:, k)' * m(:, i, ts));
                            end

                            lhs >= rhs;
                        end

                        % Increase the active time slot
                        ts = ts + 1;
                    end
            cvx_end

            m_p = m;

            % SINR
            y_p = SINR(h, m, N0);

            % Number of active time slots
            A = sum(~isinf(y_p(1, :)));

            T_m = 1 ./ (log(1 + min(y_p))); 

            R = (A + 1) / sum(T_m);

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

    [~, ~, T] = size(m);

    y = inf*ones(K, T);

    ts = 1;

    for v = nchoosek(1:K, K-1)'
        % v is the set of active UEs 

        for k = v'
            % k is the intended user 

            I = N0;

            % j goes through the interfering UEs
            for j = v'
                if (j ~= k)
                    I = I + norm(h(:, k)' * m(:, j, ts))^2;
                end
            end

            y(k, ts) = norm(h(:, k)' * m(:, k, ts))^2 / I;
        end

        % Increase the active time slot
        ts = ts + 1;
    end
end
