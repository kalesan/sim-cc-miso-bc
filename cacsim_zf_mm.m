function rate = cacsim_zf_mm(L, K, SNR_dB, varargin)
    % Zero-forcing maxmin rate.

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

        m_zf = zeros(L, msgs);

        for k = 1:msgs
            tmp = null(h(:, any(IA == k, 2))');
            [U, D, V] = svd(h'  * tmp);

            m_zf(:, k) = tmp  * V(:, 1);
            m_zf(:, k) = m_zf(:, k) / norm(m_zf(:, k));
        end

        m_p = m_zf * sqrt(SNR / msgs);

        y_p = sinr(h, m_p, N0, MA, IA);

        R_prev = inf;

        for (i = 1:ITERATIONS)
            cvx_begin
                cvx_solver('sedumi');
                cvx_quiet('true')

                variable y(K, uem) nonnegative;

                variable p(msgs, 1) nonnegative;

                variable eph;

                expression lhs(K, uem)

                maximize(eph);

                subject to
                    y >= eph;

                    % Power constraint
                    sum_square(p(:)) <= SNR;

                    % SINR Constraints
                    for k = 1:K
                        t = 1;

                        for j = 1:msgs
                            if all(MA(j, :) ~= k)
                                continue;
                            end

                            c = N0 + norm(h(:, k)' * m_p(:, j))^2;

                            lhs(k,t) = c + (c / (1 + y_p(k, t))) * (y_p(k, t) - y(k, t));

                            lhs(k,t) = lhs(k,t) - 2 * ...
                                real(m_p(:, j)' * h(:, k) * h(:, k)' * (m_p(:, j) - m_zf(:, j) * p(j)));

                            lhs(k,t) >= N0 * (1 + y_p(k, t));

                            t = t + 1;
                        end
                    end
            cvx_end

            for k = 1:msgs
                m_p(:, k) = m_zf(:, k) * p(k);
            end

            % SINR
            y_p = sinr(h, m_p, N0, MA, IA);

            % MAC rate points (K * R_MAC)
            R = [];
            for k = 1:K, for v = nchoosek(1:uem, k)'
                tmp = log(1 + sum(y_p(:, v), 2)) / k;
                R = [R tmp];
            end, end

            % Messages per UE
            M = sum(MA(:) == 1);

            R = 1 / min(R(:)); 

            rate(r, i) = (M+1) / R;

            disp(['Realization ' int2str(r) ' rate ' num2str(rate(r,i))  ...
                  ' pwr ' num2str(norm(m_p(:))^2 / SNR * 100) '%'])

            if (abs(rate(r, i) - R_prev) < 1e-8)
                rate(r, i:end) = rate(r, i);

                break;
            end

            R_prev = rate(r, i);
        end
    end
end
