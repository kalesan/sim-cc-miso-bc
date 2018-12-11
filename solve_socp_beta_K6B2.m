function m = solve_socp_beta_K6B2(h, m_p, y_p, SNR, N0, S, beta, G)
    [L, K] = size(h);

    groups = K / 2;

    K_g = K / groups;

    MAI = nchoosek(1:K_g, S);

    [~, groups, groupings] = size(m_p);

    [msgs, ~] = size(MAI);

    % Messages for each UE
    [uem] = sum(MAI(:, 1) == 1);

    m_ret = zeros(size(m_p));

    for g = 1:groupings
        cvx_begin
            cvx_solver('sedumi');

            cvx_quiet('true')

            variable y(K_g, groups) nonnegative;

            variable m(L, groups) complex;

            variable eph;

            maximize(eph);

            subject to
                % Objective in epigraph form
                y >= eph;

                % Power constraint
                sum_square_abs(m(:)) <= SNR;

                % SINR Constraints
                for f = 1:groups, for k = 1:K_g
                    h_ = h(:, G(k, g, f));
                    
                    % Interfering messages indices for k
                    IG = [1:(f-1),(f+1):groups];

                    c = N0 + sum_square_abs(vec(h_' * m_p(:, f, g)));
                    c = c + sum_square_abs(vec(h_' * m_p(:, IG, g)));

                    lhs = c / (1 + y_p(k, f, g)) ...
                          + (c / (1 + y_p(k, f, g))^2) * (y_p(k, f, g) - y(k, f));

                    % Approximation of signal power
                    lhs = lhs - (2 / (1 + y_p(k, f, g))) * ...
                        trace(real(m_p(:, f, g)' * h_ * ...
                            h_' * (m_p(:, f, g) - m(:, f))));

                    lhs = lhs - (2 / (1 + y_p(k, f, g))) * ...
                        trace(real(m_p(:, IG, g)' * h_ * ...
                            h_' * (m_p(:, IG, g) - m(:, IG))));

                    % Interference power
                    lhs >= N0 + sum_square_abs(vec(h_' * m(:, IG)));;
                end, end
        cvx_end

        m_ret(:, :, g) = m;
    end
    
    m = m_ret;
end
