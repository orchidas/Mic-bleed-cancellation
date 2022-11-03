%% test calibration algos by generating RIR with Image-source method

close all;clear all;
addpath(genpath('../room_acoustics_simulation/RIR-Generator/.'));   %rir toolbox
addpath(genpath('../BSIE_toolbox/.'));   %BSI toolbox

%% Initialization
M = 5;          % number of channels
fs = 8e3;       % sampling frequency
L = 128;        % channel length
F = 2*L;        % frame length
SNR = 40;       % signal to noise ratio (in dB)
s_seed = 1;     % seed for generating source signal
v_seed = 50;    % seed for generating noise signal
N = 5*fs;       % data length (samples)
rho = 0.8;      % step-size
lambda = 0.98;  % forgetting-factor

% Generate sensor signals and AIRs
air.c = 342;  % speed of sound
air.T60 = 0.3;  % reverberation time
air.room_dim = [5; 6; 3];  % room dimension
air.mic_spacing = 0.2;  % microphone spacing (m)
air.src_pos = [100*pi/180 2*pi/180 0.2];  % source location
air.cen_pos = [2.5; 2; 1.6];  % centre pos. of the array
[h, x, ~, s, src_pos, mic_pos, ~] = generate_data(M, L, fs, air, N, SNR, s_seed, v_seed);

% plot source-mic configuration
%figure(1);
fig = figure('Units','inches', 'Position',[0 0 4.5 2.5],'PaperPositionMode','auto');
sp = scatter3(src_pos(1), src_pos(2), src_pos(3));hold on;grid on;
sp.MarkerEdgeColor = 'b';
rp = scatter3(mic_pos(1,:), mic_pos(2,:), mic_pos(3,:), 'x');hold on;
rp.MarkerEdgeColor = 'r';
xlim([0 air.room_dim(1)]);ylim([0 air.room_dim(2)]);zlim([0 air.room_dim(3)]);
drawnow;
title('Room configuration');
legend('Source','Mics');
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');   
% print('../figures/bci_test_room.eps', '-depsc');


%%  NMCFLMS

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
% print('../figures/bci_nmcflms_cost.eps', '-depsc');

fig = figure('Units','inches', 'Position',[0 0 3.3 2.5],'PaperPositionMode','auto');
plot_filter(h(:,4),h_hat(:,4));  % filter coeff.
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');   
% print('../figures/bci_nmcflms_res.eps', '-depsc');


%% GCC-PHAT with LS

% nbins = 2048;
% nreflect = 15;
% [H,rir] = estimate_tf_gain_delay(x(:,1),x(:,4),fs,nreflect,nbins,L,h);
% 
% figure('Units','inches', 'Position',[0 0 3.3 2.5],'PaperPositionMode','auto');
% plot_filter(h(:,4), rir);
% set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');   
% % print('../figures/bci_gcc_phat.eps', '-depsc');

%% spectral ratio

frameSize = 1024;
hopSize = 512;
fftSize = frameSize*4;
win = hann(frameSize);

den = x(:,1);
DEN = get_stft_from_audio(den, frameSize, hopSize, fftSize, win);

num = x(:,4);
NUM = get_stft_from_audio(num, frameSize, hopSize, fftSize, win);

H0 =  mean(NUM./DEN,1);
rir_spec = real(ifft(H0, fftSize));

fig = figure('Units','inches', 'Position',[0 0 3.3 2.5],'PaperPositionMode','auto');
plot_filter(h(:,4), rir_spec(1:L).');
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize', 8, 'FontName','Times');   
% print('../figures/bci_spec_ratio.eps', '-depsc');

%% all methods in same plot

% time-domain fit
%
% fig = figure('Units','inches', 'Position',[0 0 3.3 2.5],'PaperPositionMode','auto');
% ValueArray = {'--';'-.';'-'};
% plt = plot(1:L, [normalize(h(:,4)), normalize(h_hat(:,4)), normalize(rir_spec(1:L)).']); grid on;
% % [plt(:).LineStyle] = ValueArray{:};
% xlabel('Time (samples)'); ylabel('Amplitude');
% legend('Simulated','NMFLMS','Spec-ratio');
% axis tight;
% set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize', 8, 'FontName','Times');   
% print('../../figures/calib_IS.eps', '-depsc');


% frequency domain fits

nfft = 2048;
freq_axis = linspace(0, fs/2, nfft/2);

fig = figure('Units','inches', 'Position',[0 0 3.3 2.5],'PaperPositionMode','auto');
semilogx(freq_axis, [get_magnitude_response(h(:,4),nfft), get_magnitude_response(h_hat(:,4),nfft), ...
    get_magnitude_response(rir_spec(1:L).',nfft)]); grid on;
xlim([20,10000]);ylim([-80,10]);
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)');
legend('Simulated','NMFLMS','Spec-ratio', 'Location', 'southeast');
axis tight;
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize', 8, 'FontName','Times');   
% print('../../figures/calib_IS_freq.eps', '-depsc');


function [h_norm] = normalize(h)
    h_norm = h./max(abs(h));
end

function [H_mag] = get_magnitude_response(h, nfft)
    H = fft(h,nfft);
    H = H(1:nfft/2);
    H = H./max(abs(H));
    H_mag = mag2db(abs(H));
end