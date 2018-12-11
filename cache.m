function [Mp, Ip] = cache(N)
    % N is the number of files


    % K is the number of UEs 
    K = N;

    % M is the memory size
    M = 1;

    Nf = reshape(1:(K*(N-M)), [N-M, K])';

    chosen = zeros(size(Nf));

    xors = zeros((K*(N-M)) / 2, 2);

    xi = 1;
    for k = 1:K
        j = k + 1;

        for l = 1:(N-M)
            if (chosen(k, l) == 1)
                continue;
            end

            ind = find(chosen(j, :) == 0);

            xors(xi, 1) = Nf(k, l);
            xors(xi, 2) = Nf(j, ind(1));
            xi = xi + 1;

            chosen(k, l) = 1;
            chosen(j, ind(1)) = 1;

            j = j + 1;
            if (k > K) 
                j = k + 1;
            end
        end
    end

    % Message allocation
    Mp = ceil(xors / (N-M));

    L = sum(Mp(:, 1) == 1);

    % Interference profile
    Ip = zeros(K, (K*(N-M) / 2) - L);

    for k = 1:K
        i = 1;
        for l = 1:(K*(N-M) / 2)
            if all(Mp(l, :) ~= k)
                Ip(k, i) = l;

                i = i + 1;
            end
        end
    end
end
