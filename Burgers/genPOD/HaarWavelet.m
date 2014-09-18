function haar = HaarWavelet(n, x0, xe, N)
    x = linspace(x0, xe, N);
    haar = ones(N, 1);
    
    if n == 1
        return
    elseif n == 2
        haar(1:floor(N/2)) = 1;
        haar(floor(N/2)+1:end) = -1;
    else
        l2 = floor(log2(n - 1));
        l3 = n - 2^l2 - 1;
        for ii = 1:N
            haar(ii) = sqrt(2^l2)*Haar_helper(2^l2*x(ii) - (l3*(xe-x0)), x0, xe);
        end
    end
    
    plot(x,haar)
end

function eval = Haar_helper(x, x0, xe)
    if x >= x0 && x <= (xe - x0) / 2
        eval = 1.0;
    elseif x > (xe - x0) / 2 && x <= xe
        eval = -1.0;
    else
        eval = 0.0;
    end
end