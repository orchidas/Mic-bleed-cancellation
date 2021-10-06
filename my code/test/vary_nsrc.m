close all;

addpath('../.');
addpath('../../kokkinis/.');
addpath(genpath('../../Eval/PEASS-Software-v2.0.1/.'));   %PEASS toolbox for evaluation

matpath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Mat files/'; 
readpath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/Num sources/';
savepath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/Num sources/';

load([matpath, 'quartet_synthesized_data_var_src.mat']);
sname = {'Va','Vcl','Vl1','Vl2'};
method = 'map'; % MAP/MLE, 
calib_type = 'spec-ratio'; %spec-ratio/gcc-phat/bci
str_add = ['_',calib_type,'_', method];  
fs = 48000;
Nsrc = [Sim(:).Nsrc];
sigma = [0];
frameSize = 1024;
hopSize = 512;
fftSize = frameSize*4;
mu = 2;

%% 

for k = 1:length(sigma)
    disp(['Testing for \sigma =' , num2str(sigma(k))]);
    if sigma(k) == 0 
        resname = ['Changing_num_sources_sigma=0',str_add,'_mu=',num2str(mu),'.mat'];
    else
        resname = ['Changing_num_sources_sigma=',num2str(sigma(k)),str_add,'.mat'];
    end

     if ~isfile([savepath,resname])

        %% MLE parameters
        calib_flag = true;     %has calibration been done?

        ml.IR_est = zeros(fftSize,2,2);
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

        
        % With known TF
        s_ideal = [];
        

        %% loop over source-mic distance

        for l = 1:length(Sim)

            disp(['Evaluating for num sources =', num2str(Sim(l).Nsrc)]);
            nmic = Sim(l).Nmic;
            nsrc = Sim(l).Nsrc;
            xmic =[];
            for i = 1:nmic
                m = Sim(l).mic(:,i);
                xmic = [xmic m];
            end  
                                    
            [s_est, H_opt] = debleed(xmic,nsrc,fs,frameSize,hopSize,fftSize,method,sigma(k),calib_flag,Sim(l).calib(1),calib_type);
            MLE(l).IR_est = H_opt;
            MLE(l).Nsrc = Sim(l).Nsrc;
            MLE(l).source = [MLE(l).source, s_est];
            
            if sigma(k) == 0     
                s_ideal = known_RTF(xmic,nsrc,Sim(l).h_ideal, frameSize, hopSize, fftSize);
                s_mwf = mcwf_dmrs(xmic, WinLen, HopSize, Wnd, beta, Kp, wdom, wresx, wresy, psdwemode);

            end
    
            %%
            for mic = 1:nmic
                
                if sigma(k) == 0
                    audiowrite([savepath, 'Ideal/', 'Nsrc=' ,num2str(Sim(l).Nsrc),'_',sname{mic},'.wav'],s_ideal(1:length(xmic),mic),fs);
%                   disp(['Playing Ideal TF separated ', sname{mic}, ' for src-mic distance ',                      ...num2str(Sim(l).dis),'m']);
%                   soundsc(s_ideal(:,mic),fs);pause(7);
                    audiowrite([savepath, 'MCWF/', 'Nsrc=' ,num2str(Sim(l).Nsrc),'_',sname{mic},'.wav'],s_mwf(1:length(xmic),mic),fs);
%                   disp(['Playing MCWF separated ', sname{mic}, ' for src-mic distance ',                          ...num2str(Sim(l).dis),'m']);
%                   soundsc(s_mwf(:,mic),fs);pause(7);
                end
                
                audiowrite([savepath, upper(method),'/sigma=',num2str(sigma(k)), '/Nsrc=' ,num2str(Sim(l).Nsrc),'_',sname{mic},str_add,'.wav'],s_est(1:length(xmic),mic),fs);

%                 disp(['Playing MLE separated ', sname{mic}, ' for src-mic distance ', num2str(Sim(l).dis),'m']);
%                 soundsc(s_est(:,mic),fs);pause(7);
              
            end

         end
            
            save([savepath, resname], 'MLE');
      


    else

        load([savepath, resname]);
        

    end


    %% evaluate objective parameters
    if sigma(k) == 0 
        resname = ['objective_results_sigma=0',str_add,'_mu=',num2str(mu),'.mat'];
    else
        resname = ['objective_results_sigma=',num2str(sigma(k)),str_add,'.mat'];
    end
    N = 4;
%     if ~isfile([savepath,resname])
        res_MLE = cell(length(Nsrc));
        res_MCWF = cell(length(Nsrc));
        res_ideal = cell(length(Nsrc));


        sdr_mle = zeros(length(Nsrc),N);
        sir_mle = zeros(length(Nsrc),N);
        sar_mle = zeros(length(Nsrc),N);
        ops_mle = zeros(length(Nsrc),N);
        tps_mle = zeros(length(Nsrc),N);
        ips_mle = zeros(length(Nsrc),N);
        aps_mle = zeros(length(Nsrc),N);

        sdr_mwf = zeros(length(Nsrc),N);
        sir_mwf = zeros(length(Nsrc),N);
        sar_mwf = zeros(length(Nsrc),N);
        ops_mwf = zeros(length(Nsrc),N);
        tps_mwf = zeros(length(Nsrc),N);
        ips_mwf = zeros(length(Nsrc),N);
        aps_mwf = zeros(length(Nsrc),N);
        
        
        sdr_ideal = zeros(length(Nsrc),N);
        sir_ideal = zeros(length(Nsrc),N);
        sar_ideal = zeros(length(Nsrc),N);
        ops_ideal = zeros(length(Nsrc),N);
        tps_ideal = zeros(length(Nsrc),N);
        ips_ideal = zeros(length(Nsrc),N);
        aps_ideal = zeros(length(Nsrc),N);

        

        for l = 1:length(Sim)

            for src = 1:Nsrc(l)
                
                originalFiles = cell(Nsrc(l),1);
                originalFiles{1} = [readpath, sname{src},'_',num2str(src),num2str(src),'_nsrc=',num2str(Sim(l).Nsrc),'.wav'];
                interf = setdiff(1:Nsrc(l), src);
                for n = 1:length(interf)
                    originalFiles{n+1} = [readpath, sname{interf(n)},'_',num2str(interf(n)),num2str(src),'_nsrc=',num2str(Sim(l).Nsrc),'.wav'];
                end
                
                estimateFileMLE = [savepath, upper(method),'/sigma=',num2str(sigma(k)), '/Nsrc=' ,num2str(Sim(l).Nsrc),'_',sname{src},str_add,'.wav'];
                options1.destDir = [savepath, upper(method),'/PEASS/'];
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

                    estimateFileMWF = [savepath, 'MCWF/', 'Nsrc=' ,num2str(Sim(l).Nsrc),'_',sname{src},'.wav'];
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

                    estimateFileIdeal = [savepath, 'Ideal/', 'Nsrc=' ,num2str(Sim(l).Nsrc),'_',sname{src},'.wav'];
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


                end

            end

        end

        save([savepath,resname],'sdr_mle','sir_mle','sar_mle','ops_mle','tps_mle','ips_mle','aps_mle',...
             'sdr_mwf','sir_mwf','sar_mwf','ops_mwf','tps_mwf','ips_mwf','aps_mwf',...
             'sdr_ideal','sir_ideal','sar_ideal','ops_ideal','tps_ideal','ips_ideal','aps_ideal');


%     end


end


%%


    
