classdef precoder_socp < precoder
    methods

    function obj = precoder_socp(varargin)
        obj = obj@precoder(varargin{:});
    end

    function rate = simulate(obj, h, S, P, N0, MAI, iterations)
        % Message allocation 
        [msgs, ~] = size(MAI);
        [uem] = sum(MAI(:, 1) == 1);

        [L, K] = size(h);

        rate = zeros(iterations, 1);

        m_p = 1/sqrt(2) * (randn(L, msgs) + randn(L, msgs)*1j);

        for k = 1:msgs
            m_p(:, k) = m_p(:, k) / norm(m_p(:, k));
        end

        m_p = m_p * sqrt(P / msgs);

        y_p = sinr(h, m_p, N0, S);

        R_prev = inf;
        
        for (i = 1:iterations)
            m = solve_socp(h, m_p, y_p, P, N0, S);
            m_p = m;

            % SINR
            y_p = sinr(h, m, N0, S);

            % MAC rate points (R_MAC)
            R = mac(y_p, uem, K);

            R = 1 / min(R(:)); 
 
            if (S == 2)
                rate(i) = (uem+1) / R;
            elseif (S == 3)
                rate(i) = (uem+3) / R;
            end

            disp(['Rate ' num2str(rate(i))  ...
                  ' P ' num2str(norm(m_p(:))^2 / P * 100) '%' ...
                  ' SNR ' num2str(10*log10(P)) 'dB'])

            if (abs(rate(i) - R_prev) < 1e-4)
                rate(i:end) = rate(i);

                break;
            end

            R_prev = rate(i);
        end
    end
    end
end
