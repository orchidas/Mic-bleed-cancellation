function [s_opt, H_opt] = mle_closed_form_batch(Nsrc, x, H_tilde, sigma)
%%
% For cost function J = \sum_{t=1}^T||x(t) - Hs(t)||^2 +
%                       1/sigma ||H_tilde - H||^2, solve for
% \nabla_s(t) J = 0, \nabla_H J = 0

%%
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
        if optimValues.fval < 1e-2  
            stop = true; 
            disp('Stopping, small enough error');
        end
    end


fun = @f; % function
nskip = size(x,1)/Nsrc;
s0 =  0.01*randn(size(x(1:nskip:end,:))) + x(1:nskip:end,:);
options = optimoptions('fsolve','Display','off', 'MaxIterations',500,...
    'MaxFunctionEvaluations',1e4, 'OutputFcn',@outfun); 
% options = optimoptions('fsolve','Display','off', 'MaxIterations',500); % optimizer options

[s_opt, fval, exitflag, output] = fsolve(fun,s0,options);
if exitflag <= 0
    disp(output.message);
end


xs_opt = x*s_opt';
ss_opt = s_opt*s_opt';
H_opt = (H_tilde + sigma*xs_opt)/(I + sigma*ss_opt);

end

