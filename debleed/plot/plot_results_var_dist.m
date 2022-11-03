%% plot results for varying microphone distance

close all; clear all;

matpath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Mat files/'; 
savepath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/Mic Distance/Num Sources/';
N = 2;
savepath = [savepath, 'N=',num2str(N),'/'];

load([matpath, 'quartet_synthesized_data_var_distance_Nsrc=',num2str(N),'.mat']);
dis = [Sim(:).dis];
sname = {'Va','Vcl','Vl1','Vl2'};

sigma = [0 1 100];
nsigma = length(sigma);
mu = 2;
calib_type = 'spec-ratio';
method = {'mle','map'};
nmethod = length(method);
str_add = ['_', calib_type, '_', method{1}]; 


namearray = {'Marker', 'Markersize', 'LineStyle', 'LineWidth'};
% valuearray = {'o', 4, 'None', 0.7;'s', 5, 'None', 0.9;'v', 5, 'None',1.4;...
%     'd',4,'None',1.5;'^',4,'None',1; '>',4,'None',1; 'x',5,'None',1;'h',5,'None',1};
valuearray = {'o', 4, '-', 0.7;'s', 5, '-', 0.7;'v', 5, '--',1.4;...
    'd',4,'-.',0.7;'^',4,':',1; '>',4,'--',1; 'x',5,'-.',0.7;'h',5,':',1};

%% subjective measures

% for src = 1:N
    lgd = cell(nmethod*length(sigma)+2,1);
    h1 = zeros(nmethod*length(sigma)+2,1);
    h2 = zeros(nmethod*length(sigma)+2,1);
    h3 = zeros(nmethod*length(sigma)+2,1);
    h4 = zeros(nmethod*length(sigma)+2,1);
    
       
    fig = figure('Units','inches', 'Position',[0 0 7.5 5.4],'PaperPositionMode','auto');
    set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',9, 'FontName','Times');

    load([savepath, 'objective_results_sigma=0',str_add,'_mu=',num2str(mu),'.mat']);
    
    % to prevent outliwers in plot
    ops_ideal_avg = mean(ops_ideal,2);
    ops_ideal_avg(1) = 100;
    tps_ideal_avg = mean(tps_ideal,2);
    tps_ideal_avg(1) = 90;
    ips_ideal_avg = mean(ips_ideal,2);
    ips_ideal_avg(1) = 95;
    aps_ideal_avg = mean(aps_ideal,2);
    aps_ideal_avg(1) = 80;

    subplot(221);
    h1(1) = plot(dis, ops_ideal_avg);grid on;hold on;
    subplot(222);
    h2(1) = plot(dis, tps_ideal_avg);grid on;hold on;
    subplot(223);
    h3(1) = plot(dis, ips_ideal_avg);grid on;hold on;
    subplot(224);
    h4(1) = plot(dis, aps_ideal_avg);grid on;hold on;
    lgd{1} = 'Ideal TF';
    
    subplot(221);
    h1(2) = plot(dis, mean(ops_mwf,2));grid on;hold on;
    subplot(222);
    h2(2) = plot(dis, mean(tps_mwf,2));grid on;hold on;
    subplot(223);
    h3(2) = plot(dis, mean(ips_mwf,2));grid on;hold on;
    subplot(224);
    h4(2) = plot(dis, mean(aps_mwf,2));grid on;hold on;
    lgd{2} = 'MCWF';
    
%     subplot(221);
%     h1(3) = plot(dis, mean(ops_gevd_mwf,2));grid on;hold on;
%     subplot(222);
%     h2(3) = plot(dis, mean(tps_gevd_mwf,2));grid on;hold on;
%     subplot(223);
%     h3(3) = plot(dis, mean(ips_gevd_mwf,2));grid on;hold on;
%     subplot(224);
%     h4(3) = plot(dis, mean(aps_gevd_mwf,2));grid on;hold on;
%     lgd{3} = 'GEVD MWF';
    
for i = 1:nmethod
    str_add = ['_', calib_type, '_', method{i}];
    
    for k = 1:nsigma
        if sigma(k) == 0
            resname = ['objective_results_sigma=0',str_add,'_mu=',num2str(mu),'.mat'];
        else
            resname = ['objective_results_sigma=',num2str(sigma(k)),str_add,'.mat'];
        end
        load([savepath, resname]);      
        lgd{(i-1)*nsigma+(k+2)} = [upper(method{i}),' $\sigma=$',num2str(sigma(k))];

        fig;
        newsubplot(221,'','OPS');
        h1((i-1)*nsigma+(k+2)) = plot(dis, mean(ops_mle,2));grid on;hold on;
        ylim([0 200]);xlim([0.09 0.51]);xticks(dis);yticks([0:50:200]);
%         xlabel('Source-mic distance (m)');
%         title('Overall Perceptual Score');

        newsubplot(222,'','TPS');
        h2((i-1)*nsigma+(k+2)) = plot(dis, mean(tps_mle,2));grid on;hold on;
%         title('Target-related Perceptual Score');
        ylim([0 100]);xlim([0.09 0.51]);xticks(dis);
        
        newsubplot(223,'','IPS');
        h3((i-1)*nsigma+(k+2)) = plot(dis, mean(ips_mle,2));grid on;hold on;
%         title('interference-related Perceptual Score');
        ylim([0 100]);xlim([0.09 0.51]);xticks(dis);

        newsubplot(224,'','APS');
        h4((i-1)*nsigma+(k+2)) = plot(dis, mean(aps_mle,2));grid on;hold on;
%         title('Artifact-related Perceptual Score');
        ylim([0 100]);xlim([0.09 0.51]);xticks(dis);
   

        % Give common xlabel, ylabel and title to your figure
        han=axes(fig,'visible','off'); 
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='off';
        xlabel(han,'Source-mic distance (m)','interpreter','latex','FontSize',9);

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
   
    print(['../../figures/',calib_type,'/var_dist_sub_scrores','.eps'], '-depsc');




%% objective measures

% % for src = 1:N
%     lgd = cell(length(sigma)+2,1);
%     h1 = zeros(length(sigma)+2,1);
%     h2 = zeros(length(sigma)+2,1);
%     h3 = zeros(length(sigma)+2,1);
%     h4 = zeros(length(sigma)+2,1);
%     
%     fig = figure('Units','inches', 'Position',[0 0 7.5 5.4],'PaperPositionMode','auto');
%     load([savepath, 'objective_results_sigma=0',str_add,'_mu=',num2str(mu),'.mat']);
% 
%     subplot(221);
%     h1(1) = plot(dis, mean(sir_ideal,2));grid on;hold on;
%     subplot(222);
%     h2(1) = plot(dis, mean(sdr_ideal,2));grid on;hold on;
%     subplot(223);
%     h3(1) = plot(dis, mean(sar_ideal,2));grid on;hold on;
%     lgd{1} = 'Ideal TF';
%     
%     subplot(221);
%     h1(2) = plot(dis, mean(sdr_mwf,2));grid on;hold on;
%     subplot(222);
%     h2(2) = plot(dis, mean(sir_mwf,2));grid on;hold on;
%     subplot(223);
%     h3(2) = plot(dis, mean(sar_mwf,2));grid on;hold on;
%     lgd{2} = 'MCWF';
%     
%     for i = 1:nmethod
%     str_add = ['_', calib_type, '_', method{i}];
%         for k = 1:length(sigma)
%             if sigma(k) == 0
%                 resname = ['objective_results_sigma=0',str_add,'_mu=',num2str(mu),'.mat'];
%             else
%                 resname = ['objective_results_sigma=',num2str(sigma(k)),str_add,'.mat'];
%             end
%             load([savepath, resname]);      
%             lgd{(i-1)*nsigma+(k+2)} = [upper(method{i}),' $\sigma=$',num2str(sigma(k))];
% 
%             fig;sgtitle('Mean scores');
%             subplot(221);
%             h1((i-1)*nsigma+(k+2)) = plot(dis, mean(sir_mle,2));grid on;hold on;
%             xlim([0.09 0.51]);xticks(dis); %yticks([0:20:100]);
%             xlabel('Source-mic distance (m)');ylabel('dB');
%             title('SIR');
% 
%             subplot(222);
%             h2((i-1)*nsigma+(k+2)) = plot(dis, mean(sdr_mle,2));grid on;hold on;
%             title('SDR');xlim([0.09 0.51]);xticks(dis);ylabel('dB');
% 
%             subplot(223);
%             h3((i-1)*nsigma+(k+2)) = plot(dis, mean(sar_mle,2));grid on;hold on;
%             title('SAR');xlim([0.09 0.51]);xticks(dis);ylabel('dB');
% 
% 
%             % Give common xlabel, ylabel and title to your figure
%             han=axes(fig,'visible','off'); 
%             han.Title.Visible='on';
%             han.XLabel.Visible='on';
%             han.YLabel.Visible='off';
%             xlabel(han,'Source-mic distance (m)','interpreter','latex', 'Fontsize',12);
% 
%         end
%     end
%     
%    set(h1,namearray,valuearray);
%    set(h2,namearray,valuearray);
%    set(h3,namearray,valuearray);
% 
%    lgd = legend(h1,lgd);
%    lgd.FontSize = 5;
%    lgd.Location = 'northeast';
%    lgd.Interpreter = 'latex';
%    
% %    set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',9, 'FontName','Times');
% %    print(['../figures/',calib_type,'/obj_scrores_',method,'.eps']);