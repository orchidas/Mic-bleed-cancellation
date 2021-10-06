close all;
clear all;

addpath('../.');
addpath('../../kokkinis/.');
addpath(genpath('../../Eval/PEASS-Software-v2.0.1/.'));   %PEASS toolbox for evaluation

matpath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Mat files/'; 
readpath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/RT60/';
savepath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/RT60/Num sources/';

N = 2;
load([matpath, 'quartet_synthesized_data_var_rt60+vol_Nsrc=',num2str(N),'.mat']);
sname = {'Va','Vcl','Vl1','Vl2'};
fs = 48000;
t60 = [Sim(:).beta];
calib_type = 'spec-ratio';  %spectral ratio/gcc-phat/BCI
method = 'map';   %MAP/MLE,
str_add = ['_', calib_type,'_', method];   
sigma = [0 1 100];
savepath = [savepath, 'N=',num2str(N),'/'];
frameSize = 1024;
hopSize = 512;
fftSize = 4*frameSize;

%% 

for k = 1:length(sigma)
    disp(['Testing for \sigma =' , num2str(sigma(k))]);
    resname = ['Changing_rt60_vol_sigma=',num2str(sigma(k)),str_add,'.mat'];

    if ~isfile([savepath,resname])

        %% MLE parameters
        calib_flag = true;     %has calibration been done?

        ml.IR_est = zeros(fftSize,N,N);
        ml.sigma = 0;
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

            disp(['Evaluating for RT60 =', num2str(Sim(l).beta), ' s']);
            xmic =[];
            for i = 1:N
                m = Sim(l).mic(:,i);
                xmic = [xmic m];
            end  
                  
            if sigma(k) == 0
                s_ideal = known_RTF(xmic, N, Sim(l).h_ideal, frameSize, hopSize, fftSize);
                s_mwf = mcwf_dmrs(xmic, WinLen, HopSize, Wnd, beta, Kp, wdom, wresx, wresy, psdwemode);
            end
            
            [s_mle, H_opt] = debleed(xmic,N,fs,frameSize,hopSize,fftSize,method,sigma(k),...
                calib_flag,Sim(l).calib(1),calib_type);
            MLE(l).IR_est = H_opt;
            MLE(l).rt60 = Sim(l).beta;
            MLE(l).V = Sim(l).V;
            MLE(l).source = [MLE(l).source, s_mle];
            

            for src = 1:N
                
                if sigma(k) == 0
                    audiowrite([savepath, 'Ideal/', 'rt60=' ,num2str(Sim(l).beta),'_',sname{src},'.wav'],s_ideal(1:length(xmic),src),fs);
                    audiowrite([savepath, 'MCWF/', 'rt60=' ,num2str(Sim(l).beta),'_',sname{src},'.wav'],s_mwf(1:length(xmic),src),fs);
                end
                audiowrite([savepath, upper(method), '/sigma=',num2str(sigma(k)), '/rt60='...
                    ,num2str(Sim(l).beta),'_',sname{src},str_add,'.wav'],...
                    s_mle(1:length(xmic),src),fs);
                
            end

        end

        save([savepath, resname], 'MLE');

    else

        load([savepath, resname]);

    end


    %% evaluate objective parameters
    resname = ['objective_results_sigma=',num2str(sigma(k)),str_add,'.mat'];

    
%     if ~isfile([savepath,resname])
        res_MLE = cell(length(t60));
        res_MCWF = cell(length(t60));
        res_ideal = cell(length(t60));
        originalFiles = cell(N);

        sdr_mle = zeros(length(t60),N);
        sir_mle = zeros(length(t60),N);
        sar_mle = zeros(length(t60),N);
        ops_mle = zeros(length(t60),N);
        tps_mle = zeros(length(t60),N);
        ips_mle = zeros(length(t60),N);
        aps_mle = zeros(length(t60),N);

        sdr_mwf = zeros(length(t60),N);
        sir_mwf = zeros(length(t60),N);
        sar_mwf = zeros(length(t60),N);
        ops_mwf = zeros(length(t60),N);
        tps_mwf = zeros(length(t60),N);
        ips_mwf = zeros(length(t60),N);
        aps_mwf = zeros(length(t60),N);
        
        sdr_ideal = zeros(length(t60),N);
        sir_ideal = zeros(length(t60),N);
        sar_ideal = zeros(length(t60),N);
        ops_ideal = zeros(length(t60),N);
        tps_ideal = zeros(length(t60),N);
        ips_ideal = zeros(length(t60),N);
        aps_ideal = zeros(length(t60),N);


        for src = 1:N

            for l = 1:length(Sim)

                originalFiles{1} = [readpath, sname{src},'_',num2str(src),num2str(src),'_rt60=',num2str(Sim(l).beta),'.wav'];
                interf = setdiff(1:N, src);
                for n = 1:length(interf)
                    originalFiles{n+1} = [readpath, sname{interf(n)},'_',num2str(interf(n)),num2str(src),'_rt60=',num2str(Sim(l).beta),'.wav'];
                end
                estimateFileMLE = [savepath, upper(method), '/sigma=',num2str(sigma(k)),...
                    '/rt60=' ,num2str(Sim(l).beta),'_',sname{src},str_add,'.wav'];
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
                    estimateFileMWF = [savepath, 'MCWF/', 'rt60=' ,num2str(Sim(l).beta),'_',sname{src},'.wav'];
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

                    estimateFileMWF = [savepath, 'Ideal/', 'rt60=' ,num2str(Sim(l).beta),'_',sname{src},'.wav'];
                    options3.destDir = [savepath, 'Ideal/PEASS/'];
                    options3.segmentationFactor = 1; % increase this integer if you experienced "out of memory" problems

                    res_ideal{l} = PEASS_ObjectiveMeasure(originalFiles,estimateFileMWF,options3);
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

end


%%

%% for the same volume, a higher T60 means more reflective surfaces (I think)
% T60 does not really affect performance of either mcwf or mle because only
% the early reflections affect the transfer function, the reverb tail is
% irrelevant. But when the volume is changed along with the T60,
    
