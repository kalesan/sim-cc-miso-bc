classdef precoder_zf < precoder
    properties
        power_loading;
    end

    methods
        function obj = precoder_zf(varargin)
            p = inputParser;
            p.KeepUnmatched = true;

            p.addOptional('power_loading', true, @islogical);

            p.parse(varargin{:});
            params = p.Results;

            obj = obj@precoder(varargin{:});

            obj.power_loading = params.power_loading;
        end

        function rate = simulate(obj, h, S, P, N0, MAI, iterations)
            [L, K] = size(h);

            [msgs, ~] = size(MAI);
            [uem] = sum(MAI(:, 1) == 1);

            rate = zeros(iterations, 1);

            m_zf = zeros(L, msgs);

            for msg = 1:msgs
                tmp = null(h(:, setdiff(1:K, MAI(msg, :)))');
                [U, D, V] = svd(h'  * tmp);

                m_zf(:, msg) = tmp  * V(:, 1);
                m_zf(:, msg) = m_zf(:, msg) / norm(m_zf(:, msg));
            end

            m_p = m_zf * sqrt(P / msgs);

            y_p = sinr(h, m_p, N0, S);

            R_prev = inf;

            for (i = 1:iterations)
                if obj.power_loading
                    m_p = solve_zf(h, m_p, m_zf, y_p, P, N0, S);
                end

                % SINR
                y_p = sinr(h, m_p, N0, S);

                % MAC rate points (K * R_MAC)
                R = [];
                for k = 1:K, for v = nchoosek(1:uem, k)'
                    tmp = log(1 + sum(y_p(:, v), 2)) / k;
                    R = [R tmp];
                end, end

                R = 1 / min(R(:)); 

                if (S == 2)
                    rate(i) = (uem+1) / R;
                elseif (S == 3)
                    rate(i) = (uem+3) / R;
                end

                disp(['Rate ' num2str(rate(i))  ...
                      ' P ' num2str(norm(m_p(:))^2 / P * 100) '%' ...
                      ' SNR ' num2str(10*log10(P)) 'dB'])

                if (abs(rate(i) - R_prev) < 1e-8)
                    rate(i:end) = rate(i);

                    break;
                end

                R_prev = rate(i);
            end
        end
    end
end
