
function [intval] = IntegrateSimpson1D(func, dp)

    N = length(func);
    intval = 0.0;

    if rem(N,2) == 1
        k_end = N/2;
    else
        k_end = (N-1)/2;
    end

    for k=1:k_end
        intval = intval + ...
            func(2*k-1) + 4.0*func(2*k) + func(2*k+1);
    end

    intval = dp/3.0*intval;
    if rem(N,2) == 0
        intval = intval + dp/2.0*(func(N-1)+func(N));
    end
end
