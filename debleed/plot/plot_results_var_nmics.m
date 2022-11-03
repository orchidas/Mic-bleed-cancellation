%% plot results for varying number of mics

close all;

matpath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Mat files/'; 
savepath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/Num mics/Num Sources/';
Nsrc = 2;
savepath = [savepath, 'N=',num2str(Nsrc),'/'];

load([matpath, 'quartet_synthesized_data_var_mics_Nsrc=',num2str(Nsrc),'.mat']);
Nmics = [Sim(:).Nmic];
sname = {'Va','Vcl','Vl1','Vl2'};

sigma = [0 1 100];
nsigma = length(sigma);
mu = 2;
calib_type = 'spec-ratio';
method = {'mle','map'};
nmethod = length(method);
str_add = ['_', calib_type, '_', method{1}]; 

namearray = {'Marker', 'MarkerSize', 'LineStyle', 'LineWidth', 'Color'};

% valuearray = {'o',4, 'None', 0.7, [0, 0.4470, 0.7410];'v', 5, 'None',1.4,[0.9290, 0.6940, 0.1250];...
%     'd',4,'None',1.5,[0.4940, 0.1840, 0.5560];'^',4,'None',1.2,[0.4660, 0.6740, 0.1880];...
%     '>',4,'None',1, [0.4660, 0.6740, 0.1880]; 'x',5,'None',1,[0.3010,  0.7450, 0.9330];...
%     'h',5,'None',1,[0.6350, 0.0780, 0.1840]};

valuearray = {'o',4, '-', 0.7, [0, 0.4470, 0.7410];'v', 5, '--',1.4,[0.9290, 0.6940, 0.1250];...
    'd',4,'-.',1,[0.4940, 0.1840, 0.5560];'^',4,':',1.2,[0.4660, 0.6740, 0.1880];...
    '>',4,'--',1, [0.3010,  0.7450, 0.9330]; 'x',5,'-.',1, [0.6350, 0.0780, 0.1840];...
    'h',5,':',1, [0, 0.4470, 0.7410]};

%% subjective scores

% for src = 1:N
    lgd = cell(nmethod*length(sigma)+1,1);
    h1 = zeros(nmethod*length(sigma)+1,1);
    h2 = zeros(nmethod*length(sigma)+1,1);
    h3 = zeros(nmethod*length(sigma)+1,1);
    h4 = zeros(nmethod*length(sigma)+1,1);
    
    load([savepath, 'objective_results_sigma=0',str_add,'_mu=',num2str(mu),'.mat']);
    
    fig = figure('Units','inches', 'Position',[0 0 7.5 5.4],'PaperPositionMode','auto');
    set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',9, 'FontName','Times');
    
    subplot(221); 
    h1(1) = plot(Nmics, mean(ips_ideal,2));grid on;hold on;
    subplot(222);h2(1) = plot(Nmics, mean(tps_ideal,2));grid on;hold on;
    subplot(223);h3(1) = plot(Nmics, mean(ips_ideal,2));grid on;hold on;
    subplot(224);h4(1) = plot(Nmics, mean(aps_ideal,2));grid on;hold on;
    lgd{1} = 'Ideal TF';
    
    
for i = 1:nmethod
      str_add = ['_', calib_type, '_', method{i}];  
    
    for k = 1:length(sigma)

        if sigma(k) == 0
            resname = ['objective_results_sigma=0', str_add,'_mu=',num2str(mu),'.mat'];
        else
            resname = ['objective_results_sigma=',num2str(sigma(k)),str_add,'.mat'];
        end
        load([savepath, resname]);      
        lgd{(i-1)*nsigma+k+1} = [upper(method{i}),' $\sigma=$',num2str(sigma(k))];

        fig;
        newsubplot(221,'','OPS'); 
        h1((i-1)*nsigma+k+1) = plot(Nmics, mean(ops_mle,2));grid on;hold on;
        ylim([0 180]);xlim([0.9, 4.1]);xticks(Nmics);
%         xlabel('Number of mics');
%         title('Overall Perceptual Score');

        newsubplot(222,'','TPS'); 
        h2((i-1)*nsigma+k+1) = plot(Nmics, mean(tps_mle,2));grid on;hold on;
%         title('Target-related Perceptual Score');
        ylim([0 100]);xlim([0.9, 4.1]);xticks(Nmics);

        newsubplot(223,'','IPS'); 
        h3((i-1)*nsigma+k+1) = plot(Nmics, mean(ips_mle,2));grid on;hold on;
%         title('interference-related Perceptual Score');
        ylim([0 100]);xlim([0.9, 4.1]);xticks(Nmics);

        newsubplot(224,'','APS');
        h4((i-1)*nsigma+k+1) = plot(Nmics, mean(aps_mle,2));grid on;hold on;
%         title('Artifact-related Perceptual Score');
        ylim([0 100]);xlim([0.9, 4.1]);xticks(Nmics);


%       Give common xlabel, ylabel and title to your figure
        han=axes(fig,'visible','off'); 
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='off';
        xlabel(han,'Number of microphones','interpreter','latex','FontSize',9);

    end
end
 
   set(h1,namearray,valuearray);
   set(h2,namearray,valuearray);
   set(h3,namearray,valuearray);
   set(h4,namearray,valuearray);

   lgd = legend(h1,lgd);
   lgd.FontSize = 5;
   lgd.Location = 'northeast';
   lgd.Interpreter = 'latex';
   
   print(['../../figures/',calib_type,'/var_mics_sub_scrores_new','.eps'], '-depsc');


%% objective scores

%     lgd = cell(length(sigma)+1,1);
%     h1 = zeros(length(sigma)+1,1);
%     h2 = zeros(length(sigma)+1,1);
%     h3 = zeros(length(sigma)+1,1);
%     h4 = zeros(length(sigma)+1,1);
%     
%     load([savepath, 'objective_results_sigma=0',str_add,'_mu=',num2str(mu),'.mat']);
%     fig = figure('Units','inches', 'Position',[0 0 3.29 2.4],'PaperPositionMode','auto');
%     subplot(221); 
%     h1(1) = plot(Nmics, mean(sir_ideal,2));grid on;hold on;
%     subplot(222);h2(1) = plot(Nmics, mean(sdr_ideal,2));grid on;hold on;
%     subplot(223);h3(1) = plot(Nmics, mean(sar_ideal,2));grid on;hold on;
%     lgd{1} = 'Ideal TF';
%     
% %     subplot(221);
% %     h1(2) = plot(Nmics, mean(sir_gevd_mwf,2));grid on;hold on;
% %     subplot(222);
% %     h2(2) = plot(Nmics, mean(sdr_gevd_mwf,2));grid on;hold on;
% %     subplot(223);
% %     h3(2) = plot(Nmics, mean(sar_gevd_mwf,2));grid on;hold on;
% %     lgd{2} = 'GEVD MWF';
%     
%     
%     for k = 1:length(sigma)
% 
%         if sigma(k) == 0
%             resname = ['objective_results_sigma=0', str_add,'_mu=',num2str(mu),'.mat'];
%         else
%             resname = ['objective_results_sigma=',num2str(sigma(k)),str_add,'.mat'];
%         end
%         load([savepath, resname]);      
% %         lgd{k+2} = ['MLE hyp =',num2str(sigma(k))];
%         lgd{k+1} = ['MAP hyp =',num2str(sigma(k))];
% 
%         subplot(221); 
%         h1(k+1) = plot(Nmics, mean(ops_mle,2));grid on;hold on;
%         xticks(Nmics);ylabel('dB');
%         title('SIR');
% 
%         subplot(222);h2(k+1) = plot(Nmics, mean(sdr_mle,2));grid on;hold on;
%         title('SDR');xticks(Nmics);ylabel('dB');
% 
%         subplot(223);h3(k+1) = plot(Nmics, mean(sar_mle,2));grid on;hold on;
%         title('SAR');xticks(Nmics);ylabel('dB');
% 
% 
% %         Give common xlabel, ylabel and title to your figure
%         han=axes(fig,'visible','off'); 
%         han.Title.Visible='on';
%         han.XLabel.Visible='on';
%         han.YLabel.Visible='off';
%         xlabel(han,'Number of microphones','interpreter','latex', 'Fontsize',14);
% 
%     end
%    set(h1,namearray,valuearray);
%    set(h2,namearray,valuearray);
%    set(h3,namearray,valuearray);
% 
%    lgd = legend(h3,lgd);
%    lgd.FontSize = 7;
%    lgd.Location = 'southeast';
   