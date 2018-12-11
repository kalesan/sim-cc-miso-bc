function rate = cacsim_socp(L, K, S, SNR_dB, varargin)
    p = inputParser;
    
    p.addOptional('realizations',  500, @isnumeric);
    p.addOptional('iterations',  100, @isnumeric);

    p.parse(varargin{:});
    params = p.Results;

    REALIZATIONS = params.realizations;
    ITERATIONS   = params.iterations;

    SNR = 10^(SNR_dB / 10);

    N0 = 1;

    % Message allocation 
    MAI = nchoosek(1:K, S);

    [msgs, ~] = size(MAI);
    [uem] = sum(MAI(:, 1) == 1);

    rate = zeros(REALIZATIONS, ITERATIONS);

    for (r = 1:REALIZATIONS)
        randn('seed', r);

        h = (1/sqrt(2)) * (randn(L, K) + 1j*randn(L, K));

        m_p = 1/sqrt(2) * (randn(L, msgs) + randn(L, msgs)*1j);

        for k = 1:msgs
            m_p(:, k) = m_p(:, k) / norm(m_p(:, k));
        end

        m_p = m_p * sqrt(SNR / msgs);

        y_p = sinr(h, m_p, N0, S);

        R_prev = inf;
        
        T = K;

        for (i = 1:ITERATIONS)
            m = solve_socp(h, m_p, y_p, SNR, N0, S);
            m_p = m;

            % SINR
            y_p = sinr(h, m, N0, S);

            % MAC rate points (R_MAC)
            R = [];
            for k = 1:K, 
                for v = nchoosek(1:uem, k)'
                    tmp = log(1 + sum(y_p(:, v), 2)) / k;
                    R = [R tmp];
                end, 
            end

            R = 1 / min(R(:)); 
 
            if (S == 2)
                rate(r, i) = (uem+1) / R;
            elseif (S == 3)
                rate(r, i) = (uem+3) / R;
            end

            disp(['Realization ' int2str(r) ' rate ' num2str(rate(r, i))  ...
                  ' pwr ' num2str(norm(m_p(:))^2 / SNR * 100) '%'])

            if (abs(rate(r, i) - R_prev) < 1e-4)
                rate(r, i:end) = rate(r, i);

                break;
            end

            R_prev = rate(r, i);
        end
    end
end
