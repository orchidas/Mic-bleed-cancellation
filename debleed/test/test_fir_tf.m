close all, clear all;

soundpath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/';  %anechoic string quartet recordings
fname = {'Va','Vcl'};
fs = 48000;
start = round(35*fs);
stop = round(42*fs);

%% loop through to create RIRs and convolve them with audio

[s1,fs] = audioread([soundpath,'Va.wav']);
s1 = s1(:,1); %convert stereo to mono
s1 = trim_audio(s1,fs,start,stop);

[s2,fs] = audioread([soundpath,'Vcl.wav']);
s2 = s2(:,1); %convert stereo to mono
s2 = trim_audio(s2,fs,start,stop);


b11 = 1;
d1 = round(0.001*fs);   %1ms delay
b12 = [1,zeros(1,d1-2),0.5];
x1 = filter(b11, 1,s1) + filter(b12,1,s2);


b22 = 1;
d2 = round(0.005*fs);   %5ms delay
b21 = [1,zeros(1,d2-2),0.3];
x2 = filter(b21,1,s1) + filter(b22,1,s2);

xmic = [x1, x2];

%% test GCC-PHAT

nbins = 2048;
imp = [1; zeros(1*fs-1,1)];
d1 = round(0.05*fs);
d2 = round(0.03*fs);
d3 = round(0.01*fs);
d4 = round(0.02*fs);
b = [zeros(1,d1-1), 0.7, zeros(1,d2-1),0.3, zeros(1,d3-1),0.1, zeros(1,d4-1),0.05].';
s1_del = filter(b,1,s1);
imp_del = filter(b,1,imp);
[alpha, tau] = estimate_tf_gain_delay(s1,s1_del,fs,3,nbins);


%%
n = 2^10;
frameSize = 1024;
hopSize = 256;
fftSize = frameSize*4;
N = 2;

h_ideal = zeros(n,N,N);
h_ideal(:,1,1) = [b11, zeros(1,n-1)];
h_ideal(:,1,2) = [b12, zeros(1,n-length(b12))];
h_ideal(:,2,1) = [b21, zeros(1,n-length(b21))];
h_ideal(:,2,2) = [b22, zeros(1,n-1)];

s_ideal = known_RTF(xmic,N,h_ideal,frameSize,hopSize,fftSize);

%%
% for src = 1:N
%     disp(['Playing separated ', fname{src}]);
%     soundsc(s_ideal(:,src),fs);pause(8);
% end


