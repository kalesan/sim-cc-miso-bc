function m = solve_socp_beta_K6B3(h, m_p, y_p, SNR, N0, S, beta)
    [L, K] = size(h);

    groups = K / 3;

    K_g = K / groups;

    MAI = nchoosek(1:K_g, S);

    groupings = nchoosek(K, K_g) / 2;

    G = zeros(K_g, groupings, groups);
    tmp = nchoosek(1:K, K_g)';
    G(:, :, 1) = tmp(:, 1:groupings);

    for g = 1:groupings
        G(:, g, 2) = setdiff(1:K, G(:,g,1));
    end

    [msgs, ~] = size(MAI);

    % Messages for each UE
    [uem] = sum(MAI(:, 1) == 1);

    cvx_begin
        cvx_solver('sedumi');

        cvx_quiet('true')

        variable y(K, uem, groups, groupings) nonnegative;

        variable m(L, msgs, groups, groupings) complex;

        variable eph(groupings, 1);

        maximize(sum(eph));

        subject to
            % Objective in epigraph form
            for g = 1:groupings, for k = 1:K_g, for v = nchoosek(1:uem, k)'
                1 + sum(y(:, v, :, g), 2) >= pow_abs(eph(g), k);
            end, end, end

            % Power constraint
            for g = 1:groupings
                sum_square_abs(vec(m(:, :, :, g))) <= SNR;
            end

            % SINR Constraints
            for g = 1:groupings, for f = 1:groups
                for k = 1:K_g
                    t = 1;

                    h_ = h(:, G(k, g, f));

                    for msg = find(any(MAI' == k))
                        % Interfering messages indices for k
                        IM = find(all(MAI' ~= k));

                        c = N0 + sum_square_abs(vec(h_' * m_p(:, [msg, IM], f, g)));

                        for f_ = [1:(f-1),(f+1):groups]
                            c = c + sum_square_abs(vec(h_' * m_p(:, :, f_, g)));
                        end

                        lhs = c / (1 + y_p(k, t, f, g)) ...
                              + (c / (1 + y_p(k, t, f, g))^2) * (y_p(k, t, f, g) - y(k, t, f, g));

                        % Approximation of signal power
                        lhs = lhs - (2 / (1 + y_p(k, t, f, g))) * ...
                            trace(real(m_p(:, [msg, IM], f, g)' * h_ * ...
                                h_' * (m_p(:, [msg, IM], f, g) - m(:, [msg, IM], f, g))));

                        for f_ = [1:(f-1),(f+1):groups]
                            lhs = lhs - (2 / (1 + y_p(k, t, f, g))) * ...
                                trace(real(m_p(:, :, f_, g)' * h_ * ...
                                    h_' * (m_p(:, :, f_, g) - m(:, :, f_, g))));
                        end

                        % Interference power
                        intf = N0 + sum_square_abs(vec(h_' * m(:, IM, f, g)));
                        for f_ = [1:(f-1),(f+1):groups]
                            intf = intf + sum_square_abs(vec(h_' * m(:, :, f_, g)));
                        end

                        lhs >= intf;

                        t = t + 1;
                    end
                end
            end, end
    cvx_end
end
