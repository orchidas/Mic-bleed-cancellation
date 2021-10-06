close all;

addpath('../.');
addpath('../../kokkinis/.');
addpath(genpath('../../Eval/PEASS-Software-v2.0.1/.'));   %PEASS toolbox for evaluation

matpath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Mat files/'; 
readpath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/Mic Distance/';
savepath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/Mic Distance/Num Sources/';

N = 2;
load([matpath, 'quartet_synthesized_data_var_distance_Nsrc=',num2str(N),'.mat']);
sname = {'Va','Vcl','Vl1','Vl2'};
fs = 48000;
calib_type = 'spec-ratio';  %spectral ratio/gcc-phat/BCI
method = 'map';   %MAP/MLE,
str_add = ['_', calib_type,'_', method];   
dis = [Sim(:).dis];
sigma = [1];
savepath = [savepath, 'N=',num2str(N),'/'];
frameSize = 1024;
hopSize = 512;
fftSize = frameSize*4;
mu = 2;

%% 

for k = 1:length(sigma)
    disp(['Testing for \sigma =' , num2str(sigma(k))]);
    if sigma(k) == 0
        resname = ['Changing_mic_distance_sigma=0',str_add,'_mu=',num2str(mu),'.mat'];
    else
        resname = ['Changing_mic_distance_sigma=',num2str(sigma(k)),str_add,'.mat'];
    end

%     if ~isfile([savepath,resname])

        %% MLE parameters
        calib_flag = true;     %has calibration been done?

        ml.IR_est = zeros(fftSize,N,N);
        ml.sigma = sigma(k);
        ml.source = [];
        MLE = repmat(ml,1,length(Sim));


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

                

        %% loop over source-mic distance

        for l = 1:length(Sim)

            disp(['Evaluating for src-mic distance =', num2str(Sim(l).dis), ' m']);
            xmic =[];
            for i = 1:N
                m = Sim(l).mic(:,i);
                xmic = [xmic m];
            end  
                        
            [s_mle, H_opt] = debleed(xmic,N,fs,frameSize,hopSize,fftSize,method,sigma(k),calib_flag,Sim(l).calib(1),calib_type);  
            MLE(l).IR_est = H_opt;
            MLE(l).dist = Sim(l).dis;
            MLE(l).source = [MLE(l).source, s_mle];
            
            if sigma(k) == 0    
                s_ideal = known_RTF(xmic, N, Sim(l).h_ideal, frameSize, hopSize, fftSize);
                s_mwf = mcwf_dmrs(xmic, WinLen, HopSize, Wnd, beta, Kp, wdom, wresx, wresy, psdwemode);
            end
    
            %%
            for src = 1:N
                
                if sigma(k) == 0
                    audiowrite([savepath, 'Ideal/', 'src_mic_dist=' ,num2str(Sim(l).dis),'_',sname{src},'.wav'],s_ideal(1:length(xmic),src),fs);
%                   disp(['Playing Ideal TF separated ', sname{mic}, ' for src-mic distance ', num2str(Sim(l).dis),'m']);
%                   soundsc(s_ideal(:,mic),fs);pause(7);
                    audiowrite([savepath, 'MCWF/', 'src_mic_dist=' ,num2str(Sim(l).dis),'_',sname{src},'.wav'],s_mwf(1:length(xmic),src),fs);
%                   disp(['Playing MCWF separated ', sname{mic}, ' for src-mic distance ', num2str(Sim(l).dis),'m']);
%                   soundsc(s_mwf(:,mic),fs);pause(7);
                end
                
                audiowrite([savepath, upper(method),'/sigma=',num2str(sigma(k)), '/src_mic_dist=' ,num2str(Sim(l).dis),'_',sname{src},str_add,'.wav'],s_mle(1:length(xmic),src),fs);
%                 disp(['Playing MLE separated ', sname{mic}, ' for src-mic distance ', num2str(Sim(l).dis),'m']);
%                 soundsc(s_mle(:,mic),fs);pause(7);
              
            end

        end
                   
        save([savepath, resname], 'MLE');


%     else
% 
%         load([savepath, resname]);
% 
%     end


    %% evaluate objective parameters
    if sigma(k) == 0
        resname = ['objective_results_sigma=0',str_add,'_mu=',num2str(mu),'.mat'];
    else
        resname = ['objective_results_sigma=',num2str(sigma(k)),str_add,'.mat'];
    end
    
%     if ~isfile([savepath,resname])
        res_MLE = cell(length(dis));
        res_MCWF = cell(length(dis));
        res_ideal = cell(length(dis));
        res_GEVD_MWF = cell(length(dis));
        originalFiles = cell(N);

        sdr_mle = zeros(length(dis),N);
        sir_mle = zeros(length(dis),N);
        sar_mle = zeros(length(dis),N);
        ops_mle = zeros(length(dis),N);
        tps_mle = zeros(length(dis),N);
        ips_mle = zeros(length(dis),N);
        aps_mle = zeros(length(dis),N);

        sdr_mwf = zeros(length(dis),N);
        sir_mwf = zeros(length(dis),N);
        sar_mwf = zeros(length(dis),N);
        ops_mwf = zeros(length(dis),N);
        tps_mwf = zeros(length(dis),N);
        ips_mwf = zeros(length(dis),N);
        aps_mwf = zeros(length(dis),N);
        
%         sdr_gevd_mwf = zeros(length(dis),N);
%         sir_gevd_mwf = zeros(length(dis),N);
%         sar_gevd_mwf = zeros(length(dis),N);
%         ops_gevd_mwf = zeros(length(dis),N);
%         tps_gevd_mwf = zeros(length(dis),N);
%         ips_gevd_mwf = zeros(length(dis),N);
%         aps_gevd_mwf = zeros(length(dis),N);
        
        sdr_ideal = zeros(length(dis),N);
        sir_ideal = zeros(length(dis),N);
        sar_ideal = zeros(length(dis),N);
        ops_ideal = zeros(length(dis),N);
        tps_ideal = zeros(length(dis),N);
        ips_ideal = zeros(length(dis),N);
        aps_ideal = zeros(length(dis),N);


        for src = 1:N

            for l = 1:length(Sim)

                originalFiles{1} = [readpath, sname{src},'_',num2str(src),num2str(src),'_mic_dist=',num2str(Sim(l).dis),'.wav'];
                interf = setdiff(1:N, src);
                for n = 1:length(interf)
                    originalFiles{n+1} = [readpath, sname{interf(n)},'_',num2str(interf(n)),num2str(src),'_mic_dist=',num2str(Sim(l).dis),'.wav'];
                end
                estimateFileMLE = [savepath, upper(method), '/sigma=',num2str(sigma(k)), '/src_mic_dist=' ,num2str(Sim(l).dis),'_',sname{src},str_add,'.wav'];
                options1.destDir = [savepath, upper(method), '/PEASS/'];
                options1.segmentationFactor = 1; % increase this integer if you experienced "out of memory" problems

                res_MLE{l} = PEASS_ObjectiveMeasure(originalFiles,estimateFileMLE,options1);
                sdr_mle(l,src) = res_MLE{l}.SDR;
                sir_mle(l,src) = res_MLE{l}.SIR;
                sar_mle(l,src) = res_MLE{l}.SAR;
                ops_mle(l,src) = res_MLE{l}.OPS;
                tps_mle(l,src) = res_MLE{l}.TPS;
                ips_mle(l,src) = res_MLE{l}.IPS;
                aps_mle(l,src) = res_MLE{l}.APS;  

                if sigma(k) == 0
                    estimateFileMWF = [savepath, 'MCWF/', 'src_mic_dist=' ,num2str(Sim(l).dis),'_',sname{src},'.wav'];
                    options2.destDir = [savepath, 'MCWF/PEASS/'];
                    options2.segmentationFactor = 1; % increase this integer if you experienced "out of memory" problems

                    res_MCWF{l} = PEASS_ObjectiveMeasure(originalFiles,estimateFileMWF,options2);
                    sdr_mwf(l,src) = res_MCWF{l}.SDR;
                    sir_mwf(l,src) = res_MCWF{l}.SIR;
                    sar_mwf(l,src) = res_MCWF{l}.SAR;
                    ops_mwf(l,src) = res_MCWF{l}.OPS;
                    tps_mwf(l,src) = res_MCWF{l}.TPS;
                    ips_mwf(l,src) = res_MCWF{l}.IPS;
                    aps_mwf(l,src) = res_MCWF{l}.APS;

                    estimateFileIdeal = [savepath, 'Ideal/', 'src_mic_dist=' ,num2str(Sim(l).dis),'_',sname{src},'.wav'];
                    options3.destDir = [savepath, 'Ideal/PEASS/'];
                    options3.segmentationFactor = 1; % increase this integer if you experienced "out of memory" problems

                    res_ideal{l} = PEASS_ObjectiveMeasure(originalFiles,estimateFileIdeal,options3);
                    sdr_ideal(l,src) = res_ideal{l}.SDR;
                    sir_ideal(l,src) = res_ideal{l}.SIR;
                    sar_ideal(l,src) = res_ideal{l}.SAR;
                    ops_ideal(l,src) = res_ideal{l}.OPS;
                    tps_ideal(l,src) = res_ideal{l}.TPS;
                    ips_ideal(l,src) = res_ideal{l}.IPS;
                    aps_ideal(l,src) = res_ideal{l}.APS;
                    
%                     estimateFileMWF = [savepath, 'GEVD-MWF/', 'src_mic_dist=' ,num2str(Sim(l).dis),'_',sname{src},'_mu=',num2str(mu),'.wav'];
%                     options2.destDir = [savepath, 'GEVD-MWF/PEASS/'];
%                     options2.segmentationFactor = 1; % increase this integer if you experienced "out of memory" problems
% 
%                     res_GEVD_MWF{l} = PEASS_ObjectiveMeasure(originalFiles,estimateFileMWF,options2);
%                     sdr_gevd_mwf(l,src) = res_GEVD_MWF{l}.SDR;
%                     sir_gevd_mwf(l,src) = res_GEVD_MWF{l}.SIR;
%                     sar_gevd_mwf(l,src) = res_GEVD_MWF{l}.SAR;
%                     ops_gevd_mwf(l,src) = res_GEVD_MWF{l}.OPS;
%                     tps_gevd_mwf(l,src) = res_GEVD_MWF{l}.TPS;
%                     ips_gevd_mwf(l,src) = res_GEVD_MWF{l}.IPS;
%                     aps_gevd_mwf(l,src) = res_GEVD_MWF{l}.APS;


                end
            end

        end

        save([savepath,resname],'sdr_mle','sir_mle','sar_mle','ops_mle','tps_mle','ips_mle','aps_mle',...
             'sdr_mwf','sir_mwf','sar_mwf','ops_mwf','tps_mwf','ips_mwf','aps_mwf',...
             'sdr_ideal','sir_ideal','sar_ideal','ops_ideal','tps_ideal','ips_ideal','aps_ideal');

%     end



end






    
