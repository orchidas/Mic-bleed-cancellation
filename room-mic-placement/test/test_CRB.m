clear all, close all, clc;

addpath('../.');
addpath('../plot/.');

fs = 8000;
room_dims = [3.5,4.6,3.2];
Nsrc = 3;
Nmic = Nsrc;
src_pos = [2.1 1.4 1.5; 1.2 1.8 2; 0.6, 0.5, 1.2];
mic_dist = 0.2;
rec_pos = src_pos - repmat([0 mic_dist 0], Nmic, 1);   %receiver positions for close mics 
rt60 = 0.5;
nsamp = 1024;
meas_noise_var = 1e-3;

% generate room impulse response matrix
[h_ir, H_ir] = generate_room_impulse_response(fs, room_dims, src_pos, rec_pos, rt60, nsamp);

%% calculate CRB for source estimation

var_bounds = calculate_CRB(H_ir, meas_noise_var);
nfft = size(H_ir,1);
fbins = linspace(0,fs/2,nfft/2);

figure;
semilogx(fbins, 20*log10(abs(var_bounds(1:nfft/2,:))));grid on;
xlabel('Frequency (Hz)'); ylabel('CRB Magnitude(dB)');