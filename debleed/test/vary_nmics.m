close all;

addpath('../.');
addpath('../../kokkinis/.');
addpath(genpath('../../Eval/PEASS-Software-v2.0.1/.'));   %PEASS toolbox for evaluation

matpath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Mat files/'; 
readpath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/Num mics/';
savepath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/Num mics/Num sources/';

Nsrc = 2;
load([matpath, 'quartet_synthesized_data_var_mics_Nsrc=',num2str(Nsrc),'.mat']);
sname = {'Va','Vcl','Vl1','Vl2'};
calib_type = 'spec-ratio';     %spectral ratio/gcc-phat
method = 'map';     % MAP/MLE, 
str_add = ['_',calib_type,'_', method];  
fs = 48000;
Nmics = [Sim(:).Nmic];
sigma = [1];
savepath = [savepath, 'N=',num2str(Nsrc),'/'];
mu = 2;

%% 

for k = 1:length(sigma)
    disp(['Testing for \sigma =' , num2str(sigma(k))]);
    if sigma(k) == 0
        resname = ['Changing_nmics_sigma=0',str_add,'_mu=',num2str(mu),'.mat'];
    else
        resname = ['Changing_nmics_sigma=',num2str(sigma(k)),str_add,'.mat'];
    end

%     if ~isfile([savepath,resname])

        %% MLE parameters
        frameSize = 1024;
        hopSize = 512;
        fftSize = frameSize*4;
        calib_flag = true;     %has calibration been done?

        ml.IR_est = zeros(fftSize,Nsrc,Nsrc);
        ml.sigma = 0;
        ml.source = [];
        MLE = repmat(ml,1,length(Sim));

        
        %% With known TF
        s_ideal = [];
        

        %% loop over Nmics

        for l = 1:length(Sim)

            disp(['Evaluating for Nmics =', num2str(Sim(l).Nmic)]);
            xmic =[];
            for i = 1:Nmics(l)*Nsrc
                m = Sim(l).mic(:,i);
                xmic = [xmic m];
            end  
                        
            if sigma(k) == 0
                s_ideal = [s_ideal,known_RTF(xmic, Nsrc, Sim(l).h_ideal, frameSize, hopSize, fftSize)];                
%                 s_gevd_mwf = debleed(xmic,Nsrc,fs,frameSize,hopSize,fftSize,'gevd_mwf',mu,calib_flag,Sim(l).calib(1), '');  
            end

            [s_mle, H_opt] = debleed(xmic,Nsrc,fs,frameSize,hopSize,fftSize,method,sigma(k),...
              'calib_flag', calib_flag, 'calib_mat', Sim(l).calib(1), 'calib_type', calib_type);
            MLE(l).IR_est = H_opt;
            MLE(l).Nmic = Sim(l).Nmic;
            MLE(l).source = [MLE(l).source, s_mle];
             

            for src = 1:Nsrc
                if sigma(k) == 0
                    audiowrite([savepath, 'Ideal/', 'nmics=' ,num2str(Sim(l).Nmic),'_',sname{src},'.wav'],s_ideal(1:length(xmic),src),fs);
%                     audiowrite([savepath, 'GEVD-MWF/', 'nmics=' ,num2str(Sim(l).Nmic),'_',sname{src},'_mu=',num2str(mu),'.wav'],s_gevd_mwf(1:length(xmic),src),fs);
                end
                audiowrite([savepath, upper(method),'/sigma=',num2str(sigma(k)), '/nmics=' ,num2str(Sim(l).Nmic),'_',sname{src},str_add,'.wav'],s_mle(1:length(xmic),src),fs);
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
        res_MLE = cell(length(Nmics));
        res_ideal = cell(length(Nmics));

        originalFiles = cell(Nsrc,1);

        sdr_mle = zeros(length(Nmics),Nsrc);
        sir_mle = zeros(length(Nmics),Nsrc);
        sar_mle = zeros(length(Nmics),Nsrc);
        ops_mle = zeros(length(Nmics),Nsrc);
        tps_mle = zeros(length(Nmics),Nsrc);
        ips_mle = zeros(length(Nmics),Nsrc);
        aps_mle = zeros(length(Nmics),Nsrc);
        
        sdr_ideal = zeros(length(Nmics),Nsrc);
        sir_ideal = zeros(length(Nmics),Nsrc);
        sar_ideal = zeros(length(Nmics),Nsrc);
        ops_ideal = zeros(length(Nmics),Nsrc);
        tps_ideal = zeros(length(Nmics),Nsrc);
        ips_ideal = zeros(length(Nmics),Nsrc);
        aps_ideal = zeros(length(Nmics),Nsrc);
       


        for src = 1:Nsrc

            for l = 1:length(Sim)

                originalFiles{1} = [readpath, sname{src},'_',num2str(src),num2str(src),'_nmics=',num2str(Sim(l).Nmic),'.wav'];
                interf = setdiff(1:Nsrc, src);
                for n = 1:length(interf)
                    originalFiles{n+1} = [readpath, sname{interf(n)},'_',num2str(interf(n)),num2str(src),'_nmics=',num2str(Sim(l).Nmic),'.wav'];
                end
                estimateFileMLE = [savepath, upper(method),'/sigma=',num2str(sigma(k)), '/nmics=' ,num2str(Sim(l).Nmic),'_',sname{src},str_add,'.wav'];
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
                    estimateFileMWF = [savepath, 'Ideal/', 'nmics=' ,num2str(Sim(l).Nmic),'_',sname{src},'.wav'];
                    options2.destDir = [savepath, 'Ideal/PEASS/'];
                    options2.segmentationFactor = 1; % increase this integer if you experienced "out of memory" problems

                    res_ideal{l} = PEASS_ObjectiveMeasure(originalFiles,estimateFileMWF,options2);
                    sdr_ideal(l,src) = res_ideal{l}.SDR;
                    sir_ideal(l,src) = res_ideal{l}.SIR;
                    sar_ideal(l,src) = res_ideal{l}.SAR;
                    ops_ideal(l,src) = res_ideal{l}.OPS;
                    tps_ideal(l,src) = res_ideal{l}.TPS;
                    ips_ideal(l,src) = res_ideal{l}.IPS;
                    aps_ideal(l,src) = res_ideal{l}.APS;

                end

            end

        end

        save([savepath,resname],'sdr_mle','sir_mle','sar_mle','ops_mle','tps_mle','ips_mle','aps_mle',...
             'sdr_ideal','sir_ideal','sar_ideal','ops_ideal','tps_ideal','ips_ideal','aps_ideal');

%     end




end




