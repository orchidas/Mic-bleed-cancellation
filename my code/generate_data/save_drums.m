close all;
clear all;

addpath(genpath('external/.'));

readpath = '../data/Drums_Noah/Take 2/';
savepath = '../data/Drums_Noah/Saved files/';
calib_path = '../data/Drums_Noah/Calibration/';
rir_path = '../data/Drums_Noah/Impulse responses/';
ext1 = '_exc.wav';
ext2 = '_1.aif';
fs = 44100;

rec_parts = {'floor tom','high hat','kick','rack tom','snare','overhead L','overhead R'};
calib_parts = {'Floor tom','Hi hat','Kick','Rack tom','Snare','Left crash','Right crash'};

nsamp = round(0.3*fs);
drums.fs = fs;
drums.Nmic = length(rec_parts);
drums.Nsrc = length(rec_parts);
drums.src_name = rec_parts;
drums.h_ideal = zeros(nsamp,drums.Nmic,drums.Nsrc);
drums.calib = struct([]);
drums.mic = [];

%% loop over mics
for nmic = 1:length(rec_parts)

    [x,fs] = audioread([readpath, rec_parts{nmic}, ext1]);
    
    %convert to mono
    x = x(:,1);
    if nmic == 1
        drums.mic = zeros(length(x),drums.Nmic);
    end
    drums.mic(:,nmic) = x;
    %normalize and save files
    audiowrite([savepath, rec_parts{nmic}, '.wav'], normalize_audio(x),fs);
end



%% setup calibration file

tau_a = 0.01;   %rise time
tau_f = 0.1;   %fall time
frac = [2,4,2,2,2,2,2]; % fraction for time constant estimation - faster decay will need smaller fraction

% for measured impulse responses
% inp_sin = audioread([rir_path,'record sine sweep_1.aif']);
% [~,sin_level] = leaky_level_detector(inp_sin,fs,tau_a,tau_a); 
% start_time = find(sin_level > 0,1);
% stop_time = length(sin_level) - find(fliplr(sin_level) > 0,1);
% inp_sin = inp_sin(start_time:stop_time);

for nsrc = 1:length(calib_parts)
    
    disp(['Calibrating ', calib_parts{nsrc}]);
    
    [x_src,fs] = audioread([calib_path,calib_parts{nsrc},'/',rec_parts{nsrc}, ext2]);
    
    %extract single drum strike
    [strikeStart, strikeEnd] = onset_detector(x_src,fs, tau_a, tau_f, frac(nsrc),calib_parts{nsrc});
    audiowrite([calib_path, calib_parts{nsrc}, '_1hit.wav'], x_src(strikeStart:strikeEnd),fs);


    for nmic = 1:length(rec_parts)
         
         [x_mic,~] =  audioread([calib_path,calib_parts{nsrc},'/',rec_parts{nmic}, ext2]);
         x_mic = x_mic(:,1);
         
         if nmic == 1
             drums.calib(1).src{nsrc}.mic = zeros(strikeEnd-strikeStart+1, drums.Nmic);
         end
         drums.calib(1).src{nsrc}.mic(:,nmic) = x_mic(strikeStart:strikeEnd);

%          disp([calib_parts{nsrc}, ' recorded by ', rec_parts{nmic},' mic.']);
%          soundsc(x_mic(strikeStart: strikeEnd),fs);pause(2);
         
            %% these measured RIRs are too noisy
            
%          [meas_sin,~] = audioread([rir_path, calib_parts{nsrc}, '/', rec_parts{nmic}, ext2]);  
%          meas_sin = meas_sin(start_time:stop_time);
%          meas_sin_cln = WienerNoiseReduction(meas_sin, fs, round(0.5*fs));   %denoise
%          [drums.h_ideal(:,nmic,nsrc)]= get_IR_from_sine_sweep(inp_sin, normalize_audio(meas_sin_cln), fs, nsamp);
%          title(['Src = ', calib_parts{nsrc}, ', Mic = ', rec_parts{nmic}]);
    end
end

save([savepath,'drums.mat'],'drums');
        
function [y_norm] = normalize_audio(y) 
    y_norm = y./max(abs(y));
end
