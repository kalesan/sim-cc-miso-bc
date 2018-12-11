function R = mac(y, uem, K)
    R = [];
    for k = 1:K, 
        for v = nchoosek(1:uem, k)'
            tmp = log(1 + sum(y(:, v), 2)) / k;
            R = [R tmp];
        end
    end
end
