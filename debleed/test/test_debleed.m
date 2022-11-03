close all;

soundpath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Mat files/';  %anechoic string quartet recordings
writepath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Results/';
sname = {'Va','Vcl'};
fs = 48000;
load([soundpath, 'quartet_synthesized_data_home_IR.mat']);

N = Sim.N;  %number of mics
%STFT parameters
frameSize = 1024;
hopSize = frameSize/2;
fftSize = frameSize*4;
win = hann(frameSize);
X = [];     %frequency-domain observation matrix of size Nframes x nbins x Nmics
calib_flag = true; %has calibration been done?
calib_type = 'gcc-phat';
method = 'mle';
xmic =[];
sigma = 1;

for i = 1:N
    m = Sim.mic(:,i);
    xmic = [xmic m];
end

s_ideal = known_RTF(xmic, N, Sim.h_direct, frameSize, hopSize, fftSize);
%%

s_mle = debleed(xmic,N,fs,frameSize,hopSize,fftSize,method,sigma, 'calib_flag', calib_flag,...
    'calib_mat', Sim, 'calib_type', calib_type);

%%
for src = 1:N    
%     disp(['Playing original ', sname{src}]);
%     soundsc(xmic(:, src), fs);pause(8);
    disp(['Playing ideal TF separated ', sname{src}]);         
    soundsc(s_ideal(:,src),fs);pause(8);
    disp(['Playing MLE separated ', sname{src}, ' for sigma =', num2str(sigma)]);
    soundsc(s_mle(:,src),fs);pause(8);

end

%% Observation - very long impulse responses introduce artifacts and perfect 
%% separation cannot be achieved with inverse filtering
