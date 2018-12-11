function m = solve_zf(h, m_p, m_zf, y_p, SNR, N0, S)
    [L, K] = size(h);

    MAI = nchoosek(1:K, S);

    [msgs, ~] = size(MAI);

    [uem] = sum(MAI(:, 1) == 1);

    cvx_begin
        cvx_solver('sedumi');
        cvx_quiet('true')

        variable y(K, uem) nonnegative;

        variable p(msgs, 1) nonnegative;

        variable eph;

        expression lhs(K, uem)

        maximize(eph);

        subject to
            % Objective in epigraph form
            for k = 1:K, for v = nchoosek(1:uem, k)'
                1 + sum(y(:, v), 2) >= pow_abs(eph, k);
            end, end

            % Power constraint
            sum_square(p(:)) <= SNR;

            % SINR Constraints
            for k = 1:K
                t = 1;

                % All messages for UE k
                for msg = find(any(MAI' == k))
                    c = N0 + norm(h(:, k)' * m_p(:, msg))^2;

                    lhs(k,t) = c + (c / (1 + y_p(k, t))) * (y_p(k, t) - y(k, t));

                    lhs(k,t) = lhs(k,t) - 2 * ...
                        real(m_p(:, msg)' * h(:, k) * h(:, k)' * (m_p(:, msg) - m_zf(:, msg) * p(msg)));

                    lhs(k,t) >= N0 * (1 + y_p(k, t));

                    t = t + 1;
                end
            end
    cvx_end

    m = zeros(size(m_p));

    for k = 1:msgs
        m(:, k) = m_zf(:, k) * p(k);
    end
end
