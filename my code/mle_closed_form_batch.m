function [s_opt, H_opt] = mle_closed_form_batch(Nsrc, x, H_tilde, sigma, varargin)
%%
% For cost function J = \sum_{t=1}^T||x(t) - Hs(t)||^2 +
%                       1/sigma ||H_tilde - H||^2, solve for
% \nabla_s(t) J = 0, \nabla_H J = 0

%%


%% parse input
p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'max_iter',500);
addParameter(p,'max_error');
addParameter(p,'max_func_count',1e6);
parse(p,varargin{:});


max_iter = p.Results.max_iter;
max_error = p.Results.max_error;
max_func_count = p.Results.max_func_count;



if sigma == 0
    H_opt = H_tilde;
    s_opt = H_opt\x;
    return;
end
    
    I = eye(Nsrc);
    % T = size(x,2);

    function y = f(s)

        %avoid loops to speed things up
        xs = x * s';
        ss = s * s';
        H = (H_tilde + sigma*xs)/(I + sigma*ss); %\nabla_H J = 0, / = A*inv(B)
        % y = (H'*H)\(H'*x) - s; % \ = inv(A)*B
        y = H\x - s; 
    end

    %optimization stopping criteria
    function stop = outfun(s,optimValues,state) 
        stop = false;
        if optimValues.fval < max_error  
            stop = true; 
            disp('Stopping, small enough error');
        end
    end


fun = @f; % function
nskip = size(x,1)/Nsrc;
s0 =  0.01*randn(size(x(1:nskip:end,:))) + x(1:nskip:end,:);


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

