function m = solve_socp(h, m_p, y_p, SNR, N0, S)
    [L, K] = size(h);

    MAI = nchoosek(1:K, S);

    [msgs, ~] = size(MAI);

    % Messages for each UE
    [uem] = sum(MAI(:, 1) == 1);

    cvx_begin
        cvx_solver('sedumi');

        cvx_quiet('true')

        variable y(K, uem) nonnegative;

        variable m(L, msgs) complex;

        variable eph;

        maximize(eph);

        subject to
            % Objective in epigraph form
            if uem == 1
                y >= eph;
            else
                for k = 1:K, for v = nchoosek(1:uem, k)'
                    1 + sum(y(:, v), 2) >= pow_abs(eph, k);
                end, end
            end

            % Power constraint
            sum_square_abs(m(:)) <= SNR;

            % SINR Constraints
            for k = 1:K
                t = 1;

                for msg = find(any(MAI' == k))
                    % Interfering messages indices for k
                    IM = find(all(MAI' ~= k));

                    c = N0 + sum_square_abs(vec(h(:, k)' * m_p(:, [msg, IM])));

                    lhs = c / (1 + y_p(k, t)) ...
                          + (c / (1 + y_p(k, t))^2) * (y_p(k, t) - y(k, t));

                    % Approximation of signal power
                    lhs = lhs - (2 / (1 + y_p(k, t))) * ...
                        trace(real(m_p(:, [msg, IM])' * h(:, k) * ...
                            h(:, k)' * (m_p(:, [msg, IM]) - m(:, [msg, IM]))));

                    % Interference power
                    lhs >= N0 + sum_square_abs(vec(h(:, k)' * m(:, IM)));;

                    t = t + 1;
                end
            end
    cvx_end
end
