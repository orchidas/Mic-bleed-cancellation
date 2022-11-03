clear all, close all, clc;

addpath('../.');
addpath('../plot/.');
addpath('../utils/.');

fs = 8000;
c = 343;
room_dims = [3.5,4.6,3.2];
src_pos = [2.1 1.4 1.5; 1.2 1.8 2; 0.5 0.6 1.3];
Nsrc = size(src_pos,1);
Nmic = Nsrc;
mic_dist = 0.2;
rec_pos = src_pos - repmat([0 mic_dist 0], Nmic, 1);   %receiver positions for close mics 
rt60 = 0.5;
nsamp = 1024;
meas_noise_var = 1e-3;

% generate room impulse response matrix
[h_ir, H_ir] = generate_room_impulse_response(fs, room_dims, src_pos, ...
    rec_pos, rt60, nsamp);
plot_room_config(room_dims, src_pos, rec_pos, 'title', 'Original mic positions');

%% estimate room mic TFs

N_extra_mics = 3;
[G, g_ir] = estimate_optimum_transfer_function(H_ir, N_extra_mics, ...
    'min_phase_flag',0);

% look at synthesized IR plots
figure;
subplot(211);
plot((0:nsamp-1)/fs, 20*log10(abs(g_ir(1:nsamp,:,2)))); hold on;
xlim([0,(nsamp-1)/fs]); ylim([-120,0]);
legend('Src 2 to room mic 1', 'Src 2 to room mic 2');
xlabel('Time');ylabel('Amplitude(dB)');
title('Optimized RIR');

% look at frequency response
col = [[0,0.447,0.741];[0.85 0.325 0.098];[0.929,0.694,0.125]];
subplot(212);
nfft = size(H_ir,1);
fbins = linspace(0,fs/2,nfft/2);
for i = 1:N_extra_mics
    semilogx(fbins, 20*log10(abs(G(1:nfft/2,i,2))), 'Color',col(i,:));
    hold on;grid on;
    semilogx(fbins, 20*log10(abs(H_ir(1:nfft/2,i,2))), 'Color',col(i,:),...
        'Linestyle','--');hold on;
end
hold off;
xlabel('Frequency (Hz)'); ylabel('CRB Magnitude(dB)');

%% look at CRBs with minimum phase reconstruction

var_bounds_close_mics = calculate_CRB(H_ir, meas_noise_var);
var_bounds_min_phase = calculate_CRB(G, meas_noise_var);

figure;
subplot(211);
semilogx(fbins, 20*log10(abs(var_bounds_close_mics(1:nfft/2,:))));
xlabel('Frequency (Hz)'); ylabel('CRB Magnitude(dB)');
xlim([20,fs/2]);
title('CRB with close mics');
subplot(212);
semilogx(fbins, 20*log10(abs(var_bounds_min_phase(1:nfft/2,:))));
xlabel('Frequency (Hz)'); ylabel('CRB Magnitude(dB)');
title('CRB with extra room mics');
xlim([20,fs/2]);


%% optimize the IRs so that they resemble real RIRs

lambda = 0.03;   %how much sparsity?
% lambda = linspace(0.03,0,nsamp).';  %sparsity decreases with time - this
% is messy
% window length for energy envelope calculation (we must be sure that h(t+beta) < h(t)
% unfortunately, this ends up determining the peak position
win_len = round(0.005 * fs);   
[h_rir,delays] = closest_viable_rir(g_ir, lambda, win_len, rt60, fs, ...
    min(room_dims), 'plot_data', 1);

% h_rir = g_ir;
% delays = find_direct_path_delay(g_ir,  min(room_dims), fs, 'plot_data', 1);


%% find optimum locations

room_mic_loc = find_optimum_mic_location(src_pos, delays, N_extra_mics, Nsrc, ...
    room_dims, mic_dist, fs);
rec_pos_new =  [rec_pos; room_mic_loc.'];
plot_room_config(room_dims, src_pos,rec_pos_new, 'title', 'Optimum mic locations');


%% generate RIRs with new locations
[h_ir_opt, H_ir_opt] = generate_room_impulse_response(fs, room_dims, src_pos, ...
    rec_pos_new, rt60, nsamp);

figure;plot([h_rir(1:nsamp/2,1,2), h_ir_opt(1:nsamp/2,1,2)]);
legend('Close mic RIR', 'Optimized RIR');

%% generate RIRs with random locations

% set seed to get same value each time
rng(123352);
rec_pos_rand = [rec_pos; 2 + 0.5*randn(N_extra_mics,3)];
plot_room_config(room_dims, src_pos, rec_pos_rand, 'marker_color','g', ...
    'title', 'Random mic locations');

[h_ir_rand, H_ir_rand] = generate_room_impulse_response(fs, room_dims, ...
    src_pos, rec_pos_rand, rt60, nsamp);

%% look at CRBs

var_bounds_close_mics = calculate_CRB(H_ir, meas_noise_var);
var_bounds_extra_mics = calculate_CRB(H_ir_opt, meas_noise_var);
var_bounds_rand_mics = calculate_CRB(H_ir_rand, meas_noise_var);


figure;
subplot(311);
semilogx(fbins, 20*log10(abs(var_bounds_close_mics(1:nfft/2,:))));
xlabel('Frequency (Hz)'); ylabel('CRB Magnitude(dB)');
xlim([20,fs/2]);
title('CRB with close mics');
subplot(312);
semilogx(fbins, 20*log10(abs(var_bounds_extra_mics(1:nfft/2,:))));
xlabel('Frequency (Hz)'); ylabel('CRB Magnitude(dB)');
title('CRB with extra room mics');
xlim([20,fs/2]);
subplot(313);
semilogx(fbins, 20*log10(abs(var_bounds_rand_mics(1:nfft/2,:))));
xlabel('Frequency (Hz)'); ylabel('CRB Magnitude(dB)');
xlim([20,fs/2]);
title('CRB with randomly placed room mics');

%% look at norms
idx = find(fbins >= 20, 1);

disp('Variance of source 1 for close mics');
norm(var_bounds_close_mics(4:end,1))
disp('Variance of source 1 for optimized room mics');
norm(var_bounds_extra_mics(4:end,1))
disp('Variance of source 1 for random room mics');
norm(var_bounds_rand_mics(4:end,1))