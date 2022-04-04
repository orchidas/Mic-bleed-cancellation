close all;
addpath('../.');
addpath('../../kokkinis/.');

readpath = '../../data/Drums_Noah/Saved files/';
load([readpath, 'drums.mat']);

src_name = drums.src_name;
Nmic = drums.Nmic;
Nsrc = drums.Nsrc;
xmic = drums.mic;
fs = 44100;
for nmic = 1:Nmic
    audiowrite([readpath, src_name{nmic}, '.wav'], xmic(:, nmic),fs);
end

sigma = 100;
method = 'map';
calib_type = 'spec_ratio';
calib_flag = true;
str_add = ['_', calib_type,'_', method];   

frameSize = 1024;
hopSize = 512;
fftSize = frameSize*4;
freqaxis = linspace(-fs/2,fs/2,fftSize);
legendcell = strings(Nsrc,Nmic);

%% MCWF parameters

WinLen = 1024;
HopSize = 512;
Wnd = hann(WinLen);

% from chapter 6, table 6.2 of his thesis
wdom = 0.6;
wresx = 0.25;
wresy = 0.15;

% I am guessing 
Kp = 5;
beta = 1;
psdwemode = 0;

        
%%
for k = 1:length(sigma)
    
%     if ~isfile([readpath,'separated_drums',str_add,'_sigma=',num2str(sigma(k)),'.mat'])
    
    % my MLE/MAP estimator
     [s_mle, H_opt, H_init] = debleed(xmic,Nsrc,fs,frameSize,hopSize,fftSize,...
         method, sigma(k),'calib_flag', calib_flag, 'calib_mat', drums.calib(1),...
         'calib_type', calib_type, 'max_error' ,1e-2, 'max_func_count' ,1e4);
     
    % Kokkinis' MCWF estimator
%     if sigma(k) == 0
%         s_mwf = mcwf_dmrs(xmic, WinLen, HopSize, Wnd, beta, Kp, wdom, wresx, wresy, psdwemode);
%     end
    
    SepDrums.IR_init = H_init;
    SepDrums.IR_est = H_opt;
    SepDrums.sigma = sigma(k);
    SepDrums.mle_source = s_mle;
    
    for nsrc = 1:Nsrc
        
        % save audio files
        audiowrite([readpath, upper(method),'/sigma=',num2str(sigma(k)),'/',...
            src_name{nsrc},str_add,'.wav'],s_mle(1:length(xmic),nsrc)./...
            abs(max(s_mle(1:length(xmic),nsrc))),fs);
%         if sigma(k) == 0
%             audiowrite([readpath, 'MCWF/',src_name{nsrc},'.wav'],...
%                 s_mwf(1:length(xmic),1,nsrc),fs);
%         end
        
%        if sigma(k) ~= 0
%             disp(['Playing recorded ', src_name{nsrc}]);
%             soundsc(xmic(:,nsrc),fs);pause(12);
%             disp(['Playing optimized ', src_name{nsrc}]);
%             soundsc(s_mle(:,nsrc),fs);pause(12);
%        end
        
    end
    
     % save results as mat file
      save([readpath,'separated_drums',str_add,'_sigma=',num2str(sigma(k)),'.mat'],'SepDrums');
%     end
end

