function [s_opt, H_opt] = map_closed_form_batch(x,H_tilde,sigma,Rss_inv,mu_s,Nsrc,sigma_w, varargin)

%%
% For cost function J = \sum_{t=1}^T||x(t) - Hs(t)||^2 + (s(t) - mu_s)^T P_s^{-1} (s(t)- mu_s) + ...
% trace((H-\tilde(H))^T.*P (H - \tilde{H})), solve for
% \nabla_s(t) J = 0, \nabla_H J = 0
%%

%% parse input
p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'max_iter', 500);
addParameter(p,'max_error',[]);
addParameter(p,'max_func_count', 1e6);
parse(p,varargin{:});


max_iter = p.Results.max_iter;
max_error = p.Results.max_error;
max_func_count = p.Results.max_func_count;


I = eye(Nsrc);
% T = size(x,2);
% Pxmu = repmat(Rss_inv*mu_s, [1,T]);

if sigma == 0
    H_opt = H_tilde;
%     s_opt = ((H_opt'*H_opt) + sigma_w*Rss_inv)\(H_opt'*x + sigma_w*Pxmu);
    s_opt = ((H_opt'*H_opt) + sigma_w*Rss_inv)\(H_opt'*x);
    return;
end


    function y = f(s)
        xs = x*s';
        ss = s*s';
        H = (H_tilde + sigma*xs)/(I + sigma*ss); %\nabla_H J = 0, / = A*inv(B)
    %     y = ((H'*H) + sigma_w*Rss_inv)\(H'*x + sigma_w*Pxmu) - s;
        y = ((H'*H) + sigma_w*Rss_inv)\(H'*x) - s;

    end

     %optimization stopping criteria
    function stop = outfun(s,optimValues,state) 
        stop = false;
        if optimValues.fval < 1e-2 
            stop = true; 
            disp('Stopping, small enough error');
        end
    end


fun = @f; % function
nskip = size(x,1)/Nsrc;
s0 = x(1:nskip:end,:) +  0.01*randn(size(x(1:nskip:end,:)));

if isempty(max_error)
    options = optimoptions('fsolve','Display','off', 'MaxIterations', max_iter, ...
          'MaxFunctionEvaluations', max_func_count); % optimizer options
else
    options = optimoptions('fsolve','Display','off', 'MaxIterations',max_iter,...
    'MaxFunctionEvaluations',max_func_count, 'OutputFcn', @outfun); 
end


[s_opt, fval, exitflag, output] = fsolve(fun,s0,options);

if exitflag <= 0
    disp(output.message);
end

xs_opt = x*s_opt';
ss_opt = s_opt*s_opt';
H_opt = (H_tilde + sigma*xs_opt)/(I + sigma*ss_opt);


end


