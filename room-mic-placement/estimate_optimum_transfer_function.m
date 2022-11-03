function [G, g_ir] = estimate_optimum_transfer_function(H_ir, N_extra_mics, varargin)

% Optimum transfer function for minimum variance with extra mic
% placement.

p = inputParser;
p.KeepUnmatched = true;
% addParameter(p,'optimize_flag',0);
addParameter(p, 'min_phase_flag',0);
parse(p,varargin{:});
% optimize_flag = p.Results.optimize_flag;
min_phase = p.Results.min_phase_flag;


[Nfreq,Nmics,Nsrc] = size(H_ir);
if (N_extra_mics > Nmics)
    warning('Number of extra mics cannot be greater than sources!');
    N_extra_mics = Nmics;
end
G = zeros(Nfreq, N_extra_mics, Nsrc);

for k = 1:Nfreq
    H_freq = reshape(H_ir(k,:,:), [Nmics, Nsrc]);
    [U,D,V] = eig(H_freq'*H_freq);
    if (N_extra_mics < Nmics)
        [vals, idx] = sort(abs(diag(D)));
        D = diag(vals(1:N_extra_mics));
        U = U(:,idx(1:N_extra_mics));
    end
    G(k,:,:) = sqrt(eye(N_extra_mics) - D)*U';

end

% if optimize_flag
%     %% estimate modal frequencies
% 
%     nsrc = Nsrc-1;
%     win_len = 1;
%     mode_idx = estimate_mode_locations(H_ir,nsrc,win_len);
% 
%     %% re-do optimization with inequality constraints if the bin 
%     % corresponds to a mode frequency
% 
%     G_opt=zeros(size(G));
% 
%     for k = 1:Nfreq/2
%         %some constraints will need to be added 
%         if ismember(k, mode_idx)
%             H_cur = reshape(H_ir(k,:,:), [Nmics,Nsrc]);
%             G_prev = reshape(G(k-1,:,:), [N_extra_mics, Nsrc]);
%             G_next = reshape(G(k+1,:,:),[N_extra_mics, Nsrc]);
%             G_cur =  reshape(G(k,:,:),[N_extra_mics, Nsrc]);

%             %cvx solver
%             G_opt(k,:,:) = cvx_solve_matrix(abs(G_prev), abs(G_next), G_cur);
%         else
%             G_opt(k,:,:) = G(k,:,:);
%         end
%     end
%     G_opt(Nfreq/2+1:end,:,:) = fliplr(conj(G_opt(1:Nfreq/2,:,:))); %conjugate symmetric spectrum 
% else
%     G_opt = G;
% end

%% inverse FFT to get time domain response

if min_phase
    g_ir = minimum_phase_IR(abs(G));
    G = fft(g_ir,Nfreq,1);

else
    g_ir = ifft(G, Nfreq, 1, 'symmetric');
end


end