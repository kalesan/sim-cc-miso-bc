function y = sinr(h, m, N0, S)
    [L, K] = size(h);
    [~, msgs] = size(m);

    MAI = nchoosek(1:K, S);

    y = zeros(K, sum(MAI(:, 1) == 1));

    for k = 1:K
        t = 1;

        for msg = find(any(MAI' == k))
            y(k, t) = norm(h(:, k)' * m(:, msg))^2 / ...
                        (N0 + norm(h(:, k)' * m(:, find(all(MAI' ~= k))))^2);

            t = t + 1;
        end
    end
end
