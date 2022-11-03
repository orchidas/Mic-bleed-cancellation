function [H,rir] = estimate_tf_gain_delay(x1,x2,fs,nreflect,nbins,L,h_meas)
%
% Approximate transfer function between two mics as a gain and delay term
% (can also include early reflections - nreflect)
% Inputs
% x1,x2 - signals from 2 mics
% nreflect - number of early reflections to detect
% nbins - number of frequency bins in estimated transfer function
% L - length of RIR (optional)
%%

if nargin < 4
    nreflect = 0;
end

maxpeaks = 1+nreflect;  %direct path + early reflections
tau = gcc_phat(x1,x2,fs,maxpeaks);  %time delays estimated with GCC-PHAT

%TF estimate
l1 = length(x1);
l2 = length(x2);
nfft = 2^nextpow2(l1+l2);
h = fft(x2, nfft)./fft(x1,nfft);
h = h(1:nfft/2);
omega = linspace(0,pi,nfft/2).';

%matrix formed with time delays
M = exp(-1i*(omega*(tau*fs)'));
alpha = real(M\h);    %gains

omega_nbins = linspace(-pi,pi,nbins).';
M_nbins = exp(-1i*(omega_nbins*(tau*fs)'));
H = M_nbins*alpha;
H = fftshift(H);

if nargin > 6
    rir = zeros(L,1);
    rir(round(tau*fs)+1) = alpha;
%     figure; plot_filter(h_meas,rir);xlim([0, 500]);  % filter coeff.
%     the delay between original and estimated is because we are
%     approximating the source with what is captured by the closest mic.In
%     reality, there is always a finite delay between those

end

end

