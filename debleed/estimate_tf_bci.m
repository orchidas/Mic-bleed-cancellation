function [h_hat, H] = estimate_tf_bci(x, fs, M, L, nbins, h)
%% estimate transfer function with BCI algorithms
% Inputs
% x - multichannel mic signal
% fs - sampling rate
% M - number of mics
% L - length of channel (RIR)
% nbins - number of frequency bins
% h - actual RIR
%%

addpath(genpath('../BSIE_toolbox/.'));   %BSI toolbox

F = 2*L;        % frame length
N = length(x);       % data length (samples)
rho = 0.2;      % step-size
lambda = 0.98;  % forgetting-factor

% Initialize NMCFLMS
%[h_hat, P_k_avg, Pn] = init_nmcflms_sc(L, F, M, x(1:F,:));
[h_hat, P_k_avg] = init_rnmcflms(L, F, M, x(1:F,:));
% ns = F-L+1;
ns = fix(L/8);
B = fix(N/ns);  % number of input blocks of length L
J = zeros(N,1);
npm_dB = zeros(B,1);




%% Processing Loop: run NMCFLMS [2]

%wbar = waitbar(0,'NMCFLMS');
for bb = 1 : B
    %waitbar(bb/B);
    if bb == 1
        xm = [zeros(F-ns,M); x(1:ns,:)]; 
    else
        xm = [xm(ns+1:end,:); x((bb-1)*ns+1:bb*ns,:)]; 
    end
    Xm = fft(x(1:F,:),F);
    P_x = conj(Xm).*Xm;
    delta = (M-1)*mean(mean(P_x));

    %[h_hat, P_k_avg, Pn, J_tmp] = nmcflms_sc(xm, h_hat,...
    %    P_k_avg, Pn, rho, lambda, delta);
    %    J((bb-1)*ns+1:(bb)*ns) = J_tmp(F-L-ns+2:end);

    [h_hat, P_k_avg] = rnmcflms(xm, h_hat, P_k_avg,...
        rho, lambda, delta);
    
    if nargin > 5
        npm_dB(bb) = 20*log10(npm(h, h_hat));  
    end
end
%close(wbar);
H = fft(h_hat, nbins);

%% plot results

% figure(1); plot_J(J, fs);  % cost function
% legend('NMCFLMS');
% title(['L= ',num2str(L), ', \rho= ',num2str(rho), ...
%     ', \lambda= ',num2str(lambda)]);
% 
% if nargin > 5
%     figure(2); plot_npm(npm_dB, fs, ns);  % NPM
%     legend('NMCFLMS');
%     title(['L= ',num2str(L), ', \rho= ',num2str(rho), ...
%         ', \lambda= ',num2str(lambda)]);
% 
%     figure(3); plot_filter(h,h_hat);  % filter coeff.
% end

end

