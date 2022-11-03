%% plot results for varying number of sources

close all;

matpath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Mat files/'; 
savepath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/Num sources/';
figsavepath = 'figures/';

load([matpath, 'quartet_synthesized_data_var_src.mat']);
Nsrc = [Sim(:).Nsrc];
sname = {'Va','Vcl','Vl1','Vl2'};

sigma = [0 1 100];
nsigma = length(sigma);
mu = 2;
calib_type = 'spec-ratio';
method = {'mle','map'};
nmethod = length(method);
str_add = ['_', calib_type, '_', method{i}]; 


namearray = {'Marker', 'Markersize', 'LineStyle', 'LineWidth'};
% valuearray = {'o', 4, 'None', 0.7;'s', 5, 'None', 0.9;'v', 5, 'None',1.4;...
%     'd',4,'None',1.5;'^',4,'None',1; '>',4,'None',1; 'x',5,'None',1;'h',5,'None',1};
valuearray = {'o', 4, '-', 0.7;'s', 5, '-', 0.7;'v', 5, '--',1.4;...
      'd',4,'-.',1;'^',4,':',1; '>',4,'--',1; 'x',5,'-.',1;'h',5,':',1};


%% subjective measures

% for src = 1:N
    lgd = cell(nmethod*length(sigma)+2,1);
    h1 = zeros(nmethod*length(sigma)+2,1);
    h2 = zeros(nmethod*length(sigma)+2,1);
    h3 = zeros(nmethod*length(sigma)+2,1);
    h4 = zeros(nmethod*length(sigma)+2,1);
    
    load([savepath, ['objective_results_sigma=0',str_add,'_mu=',num2str(mu),'.mat']]);

    fig = figure('Units','inches', 'Position',[0 0 7.5 5.4],'PaperPositionMode','auto');
    set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',9, 'FontName','Times');
    
    subplot(221);
    h1(1) = plot(Nsrc, sum(ops_ideal,2)./Nsrc');grid on;hold on;
    subplot(222);
    h2(1) = plot(Nsrc, sum(tps_ideal,2)./Nsrc');grid on;hold on;
    subplot(223);
    h3(1) = plot(Nsrc, sum(ips_ideal,2)./Nsrc');grid on;hold on;
    subplot(224);
    h4(1) = plot(Nsrc, sum(aps_ideal,2)./Nsrc');grid on;hold on;
    lgd{1} = 'Ideal TF';
    
    subplot(221);
    h1(2) = plot(Nsrc, sum(ops_mwf,2)./Nsrc');grid on;hold on;
    subplot(222);
    h2(2) = plot(Nsrc, sum(tps_mwf,2)./Nsrc');grid on;hold on;
    subplot(223);
    h3(2) = plot(Nsrc, sum(ips_mwf,2)./Nsrc');grid on;hold on;
    subplot(224);
    h4(2) = plot(Nsrc, sum(aps_mwf,2)./Nsrc');grid on;hold on;
    lgd{2} = 'MCWF';
    
for i = 1:nmethod
    str_add = ['_', calib_type, '_', method{i}];
    
    for k = 1:nsigma
        if sigma(k) == 0
            resname = [ 'objective_results_sigma=0',str_add,'_mu=',num2str(mu),'.mat'];
        else
            resname = ['objective_results_sigma=',num2str(sigma(k)),str_add,'.mat'];
        end
        load([savepath, resname]);      
        lgd{(i-1)*nsigma+(k+2)} = [upper(method{i}),' $\sigma=$',num2str(sigma(k))];

        fig;
        newsubplot(221,'','OPS');
        h1((i-1)*nsigma+(k+2)) = plot(Nsrc, sum(ops_mle,2)./Nsrc');grid on;hold on;
        ylim([0 150]);xlim([1.9,4.1]);xticks(Nsrc);yticks([0:50:200]);
%         ylabel('APS');  xlabel('Number of sources');
%          title('Overall Perceptual Score');

        newsubplot(222,'','TPS');
        h2((i-1)*nsigma+(k+2)) = plot(Nsrc, sum(tps_mle,2)./Nsrc');grid on;hold on;
%         title('Target-related Perceptual Score');
        ylim([0 100]);xlim([1.9,4.1]);xticks(Nsrc);
        
        newsubplot(223,'','IPS');
        h3((i-1)*nsigma+(k+2)) = plot(Nsrc, sum(ips_mle,2)./Nsrc');grid on;hold on;
%         title('Interference-related Perceptual Score');
        ylim([0 100]);xlim([1.9,4.1]);xticks(Nsrc);

        newsubplot(224,'','APS');
        h4((i-1)*nsigma+(k+2)) = plot(Nsrc, sum(aps_mle,2)./Nsrc');grid on;hold on;
%         title('Artifact-related Perceptual Score');
        ylim([0 100]);xlim([1.9,4.1]);xticks(Nsrc);


        % Give common xlabel, ylabel and title to your figure
        han=axes(fig,'visible','off'); 
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='off';
        xlabel(han,'Number of sources','interpreter','latex', 'Fontsize',9);

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

   
print(['../../figures/',calib_type,'/var_src_sub_scrores','.eps'], '-depsc');
%  set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',9, 'FontName','Times');
%  print('../figures/var_nsrc_N=2_Aps.eps', '-depsc');



%% objective measures

% lgd = cell(length(sigma)+2,1);
% h1 = zeros(length(sigma)+2,1);
% h2 = zeros(length(sigma)+2,1);
% h3 = zeros(length(sigma)+2,1);
% h4 = zeros(length(sigma)+2,1);
% 
% load([savepath, ['objective_results_sigma=0',str_add,'_mu=',num2str(mu),'.mat']]);
% fig = figure;
% %     fig = figure('Units','inches', 'Position',[0 0 3.29 2.4],'PaperPositionMode','auto');
% subplot(221);
% h1(1) = plot(Nsrc, sum(sir_ideal,2)./Nsrc');grid on;hold on;
% subplot(222);
% h2(1) = plot(Nsrc, sum(sdr_ideal,2)./Nsrc');grid on;hold on;
% subplot(223);
% h3(1) = plot(Nsrc, sum(sar_ideal,2)./Nsrc');grid on;hold on;
% lgd{1} = 'Ideal TF';
% 
% subplot(221);
% h1(2) = plot(Nsrc, sum(sir_mwf,2)./Nsrc');grid on;hold on;
% subplot(222);
% h2(2) = plot(Nsrc, sum(sdr_mwf,2)./Nsrc');grid on;hold on;
% subplot(223);
% h3(2) = plot(Nsrc, sum(sar_mwf,2)./Nsrc');grid on;hold on;
% lgd{2} = 'MCWF';
% 
% 
% for k = 1:length(sigma)
%     if sigma(k) == 0
%         resname = [ 'objective_results_sigma=0',str_add,'_mu=',num2str(mu),'.mat'];
%     else
%         resname = ['objective_results_sigma=',num2str(sigma(k)),str_add,'.mat'];
%     end
%     load([savepath, resname]);      
%      lgd{k+2} = ['MLE hyp =',num2str(sigma(k))];
% %         lgd{k+3} = ['MAP hyp =',num2str(sigma(k))];
% 
%     fig;sgtitle('Mean scores');
%     subplot(221);
%     h1(k+2) = plot(Nsrc, sum(sir_mle,2)./Nsrc');grid on;hold on;
%     xticks(Nsrc);ylabel('dB');
% %         xlabel('Number of sources');
%     title('SIR');
% 
%     subplot(222);
%     h2(k+2) = plot(Nsrc, sum(sdr_mle,2)./Nsrc');grid on;hold on;
%     title('SDR');xticks(Nsrc);ylabel('dB');
% 
%     subplot(223);
%     h3(k+2) = plot(Nsrc, sum(sar_mle,2)./Nsrc');grid on;hold on;
%     title('SAR');xticks(Nsrc);ylabel('dB');
% 
% 
%     % Give common xlabel, ylabel and title to your figure
%     han=axes(fig,'visible','off'); 
%     han.Title.Visible='on';
%     han.XLabel.Visible='on';
%     han.YLabel.Visible='off';
%     xlabel(han,'Number of sources','interpreter','latex', 'Fontsize',9);
% 
% end
% set(h1,namearray,valuearray);
% set(h2,namearray,valuearray);
% set(h3,namearray,valuearray);
% 
% lgd = legend(h1,lgd);
% lgd.FontSize = 7;
% lgd.Location = 'northeast';
