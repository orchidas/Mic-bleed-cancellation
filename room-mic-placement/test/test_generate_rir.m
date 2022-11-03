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
nsamp = 512;

plot_room_config(room_dims, src_pos, rec_pos);
h_ir = generate_room_impulse_response(fs, room_dims, src_pos, rec_pos, rt60, nsamp);

%%
figure;
plot((0:nsamp-1)/fs, 20*log10(abs(h_ir(:,:,2))));
xlim([0,(nsamp-1)/fs]); ylim([-120,0]);
legend('Src 2 to mic 1', 'Src 2 to mic 2', 'Src 2 to mic 3');
xlabel('Time');ylabel('Amplitude(dB)');
title('Synthesized RIR');
