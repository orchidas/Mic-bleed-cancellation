function [h] = cvx_optimize_rir(g, lambda, K, time_const, fs)
%% optimizie the solution obtained from min variance objective
% so that it resembles an RIR
% min ||h - g||^2 + \lambda*|h| (for sparsity) subject to
% h(n+k) = exp(-K/fs*time_const) * h(n)

% the increase in RIR amplitude with time is caused by the fact that I 
% cannot enforce A*abs(h) < 0 in DCP
N = length(g);
e = ones(N-1,1);
A = spdiags([exp(-K/(fs*time_const)).*e -e], [0 K-1], N-K, N);

cvx_begin
    variable h(N);
    minimize(norm((h - abs(g))) + norm(lambda.*h,1));
    subject to
        A*h == zeros(N-K,1);   %decreasing RIR
cvx_end

% h0 = mean(g) + var(g)*randn(N,1);
% [h, fval, exitflag, output] = fmincon(@costfn, h0, -A, zeros(N-tau,1),[],[], min(g), max(g));
% 
%     function err = costfn(h)
%         err = norm((g - h)) + lambda*norm(h,1);
%     end

end

