function [var_bounds] = calculate_CRB(H, meas_noise_var)
%% 
% find fisher information matrix and cramer rao bounds from transfer function matrix
% at each frequency bin

[Nfreq, Nmic, Nsrc] =  size(H);
var_bounds = zeros(Nfreq, Nsrc);
var_bounds_est = zeros(Nfreq, Nsrc);

for k = 1:Nfreq
    H_freq = reshape(H(k,:,:), [Nmic, Nsrc]);
    fisher_mat = (H_freq'*H_freq);
    var_bounds(k,:) = meas_noise_var*diag(inv(fisher_mat));
    
    %alternate calculation of bounds, these should match
    for i = 1:Nsrc
        % principal angle between each column and matrix without that column
        H_comp = H_freq;
        H_comp(:,i) = [];
    	psi = subspace(H_freq(:,i), H_comp);
        var_bounds_est(k,i) = meas_noise_var/(norm(H_freq(:,i))^2 * sin(psi)^2);
    end
    
    err = sum((20*log10(var_bounds(k,:)) - 20*log10(var_bounds_est(k,:))).^2)/Nsrc;
    assert(err < 1e-6);

end


