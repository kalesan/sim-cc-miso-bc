function rate = cacsim_sinr(L, K, SNR_dB, varargin)
    p = inputParser;
    
    p.addOptional('realizations',  500, @isnumeric);
    p.addOptional('iterations',  100, @isnumeric);

    p.parse(varargin{:});
    params = p.Results;

    REALIZATIONS = params.realizations;
    ITERATIONS   = params.iterations;

    SNR = 10^(SNR_dB / 10);

    N0 = 1;

    T = K;

    [MA, IA] = cache(K);

    [msgs, ~] = size(MA);
    [uem] = sum(MA(:, 1) == 1);

    rate = zeros(REALIZATIONS, ITERATIONS);

    for (r = 1:REALIZATIONS)
        randn('seed', r);

        h = (1/sqrt(2)) * (randn(L, K) + 1j*randn(L, K));

        m_p = 1/sqrt(2) * (randn(L, msgs) + randn(L, msgs)*1j);

        m_p = m_p / (norm(m_p(:))) * sqrt(SNR);

        y_p = SINR(h, m_p, N0, MA, IA);

        R_prev = inf;
        
        T = K;

        for (i = 1:ITERATIONS)
            cvx_begin
                cvx_solver('sedumi');

                cvx_quiet('true')

                variable y(msgs, K) nonnegative;

                variable m(L, msgs) complex;

                variable eph(msgs, 1);

                maximize(sum(eph));

                subject to
                    % Objective in epigraph form
                    for msg = 1:msgs
                        y(msg, :) >= eph(msg);
                    end

                    % Power constraint
                    for msg = 1:msgs
                        sum_square_abs(m(:, msg)) <= SNR;
                    end

                    % SINR Constraints
                    for msg = 1:msgs
                        for k = MA(msg, :)
                            c = N0 + norm(h(:, k)' * m_p(:, msg))^2;

                            lhs = c / (1 + y_p(msg, k)) ...
                                  + (c / (1 + y_p(msg, k))^2) * (y_p(msg, k) - y(msg, k));

                            lhs = lhs - (2 / (1 + y_p(msg, k))) * ...
                                real(m_p(:, msg)' * h(:, k) * h(:, k)' * (m_p(:, msg) - m(:, msg)));

                            lhs >= N0;
                        end
                    end, 
            cvx_end

            m_p = m;

            % SINR
            y_p = SINR(h, m, N0, MA, IA);

            M = sum(MA(:) == 1);

            % Rate for each message
            T_m = 1 ./ log(1 + min(y_p')); 

            R = (M + 1) / sum(T_m);

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

function y = SINR(h, m, N0, MA, IA)
    [L, K] = size(h);
    [~, msgs, T] = size(m);

    [uem] = sum(MA(:, 1) == 1);
    y = inf*ones(msgs, K);

    for msg = 1:msgs;
        for k = MA(msg, :)
            y(msg, k) = norm(h(:, k)' * m(:, msg))^2 / N0;
        end
    end
end
