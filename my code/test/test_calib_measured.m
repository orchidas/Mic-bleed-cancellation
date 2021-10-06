%% test calibration methods with recorded RIRs
addpath(genpath('../../BSIE_toolbox/.'));   %BSI toolbox
close all;

%% load recorded RIRs
soundpath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Mat files/';  %anechoic string quartet recordings
load([soundpath, 'quartet_synthesized_data_home_IR.mat']);
M = Sim.N;  %number of sources
fs = 48000;
N = 5*fs;       % data length (samples)
SNR = 40;       % signal to noise ratio (in dB)
s_seed = 1;     % seed for generating source signal
v_seed = 50;    % seed for generating noise signal
L = 128;         %length of RIR
F = 2*L;        % frame length
rho = 0.8;      % step-size
lambda = 0.98;  % forgetting-factor


h = reshape(Sim.h_ideal(1:L,:,1), [L,M]);   %TF from src 1 to all receivers
figure; plot(h(:,2));xlim([0,500]);ylim([-0.5,0.5]);
grid on;

rng(s_seed);
s = randn(N, 1);    %generate white noise
z = zeros(N,M);    %filter with impulse response
for i = 1:M
    z(:,i) = fftfilt(h(:,i),s);
end
rng(v_seed);
v = randn(N,M);
v = sqrt(var(s)*norm(h(:),2).^2 / (10^(SNR/10)*M*mean(var(v)))).*v;
x = z + v;      %noisy sensor signals

%% NMCFLMS

% Initialize
[h_hat, P_k_avg, Pn] = init_nmcflms_sc(L, F, M, x(1:F,:));
ns = F-L+1;
B = fix(N/ns);  % number of input blocks of length L
J = zeros(N,1);
npm_dB = zeros(B,1);

Xm = fft(x(1:F,:),F);
P_x = conj(Xm).*Xm;
delta = (M-1)*mean(mean(P_x));

% Processing Loop: run NMCFLMS [2]
wbar = waitbar(0,'NMCFLMS');
for bb = 1 : B
    waitbar(bb/B);
    if bb == 1
        xm = [zeros(F-ns,M); x(1:ns,:)]; 
    else
        xm = [xm(ns+1:end,:); x((bb-1)*ns+1:bb*ns,:)]; 
    end
    [h_hat, P_k_avg, Pn, J_tmp] = nmcflms_sc(xm, h_hat,...
        P_k_avg, Pn, rho, lambda, delta);
    J((bb-1)*ns+1:(bb)*ns) = J_tmp(F-L-ns+2:end);
    npm_dB(bb) = 20*log10(npm(h, h_hat));   
end
close(wbar);

% plot results
% figure(2); plot_npm(npm_dB, fs, ns);  % NPM
% legend('NMCFLMS');
% title(['L= ',num2str(L), ', \rho= ',num2str(rho), ...
%     ', \lambda= ',num2str(lambda), ', SNR= ', ...
%     num2str(SNR), ', M= ', num2str(M)]);

fig = figure('Units','inches', 'Position',[0 0 3.3 2.5],'PaperPositionMode','auto');
plot_J(J, fs);  % cost function
legend('NMCFLMS');
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');   
% print('../../figures/bci/bci_calib_meas_cost.eps', '-depsc');

fig = figure('Units','inches', 'Position',[0 0 3.3 2.5],'PaperPositionMode','auto');
plot_filter(h(:,2),h_hat(:,2));  % filter coeff.
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');   
% print('../../figures/bci/bci_calib_meas_res.eps', '-depsc');

%% GCC-PHAT

nbins = 2048;
nreflect = 50;
[H,rir] = estimate_tf_gain_delay(x(:,1),x(:,2),fs,nreflect,nbins,L,h);

figure('Units','inches', 'Position',[0 0 3.3 2.5],'PaperPositionMode','auto');
plot_filter(h(:,2), rir);xlim([1,L]);
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');   
print('../../figures/gcc-phat/gcc-phat_calib_meas.eps', '-depsc');

%% spectral ratio

frameSize = 1024;
hopSize = 512;
fftSize = frameSize*4;
win = hann(frameSize);

den = x(:,1);
DEN = get_stft_from_audio(den, frameSize, hopSize, fftSize, win);

num = x(:,2);
NUM = get_stft_from_audio(num, frameSize, hopSize, fftSize, win);

H0 =  mean(NUM./DEN,1);
rir_spec = real(ifft(H0, fftSize));

fig = figure('Units','inches', 'Position',[0 0 3.3 2.5],'PaperPositionMode','auto');
plot_filter(h(:,2), rir_spec(1:L).');
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize', 8, 'FontName','Times');   
% print('../../figures/spec-ratio/spec-ratio_calib_meas.eps', '-depsc');


