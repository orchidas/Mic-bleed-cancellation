%% plot results for varying rt60 and volume


close all;clear all;

matpath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Mat files/'; 
savepath = '../../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/RT60/Num Sources/';
N = 2;
savepath = [savepath, 'N=',num2str(N),'/'];

load([matpath, 'quartet_synthesized_data_var_rt60+vol_Nsrc=',num2str(N),'.mat']);
rt60 = [Sim(:).beta];
sname = {'Va','Vcl','Vl1','Vl2'};

sigma = [0 1 100];
nsigma = length(sigma);
mu = 2;
calib_type = 'spec-ratio';
method = {'mle', 'map'};
nmethod = length(method);
str_add = ['_', calib_type, '_', method{1}]; 


namearray = {'Marker', 'Markersize', 'LineStyle', 'LineWidth'};
% valuearray = {'o', 4, 'None', 0.7;'s', 5, 'None', 0.9;'x', 5, 'None',1.4;'d',4,'none',1.5;'^',4,'None',1;}; %'>','-.',1};
valuearray = {'o', 4, 'None', 0.7;'s', 5, 'None', 0.9;'v', 5, 'None',1.4;...
    'd',4,'None',1.5;'^',4,'None',1; '>',4,'None',1; 'x',5,'None',1;'h',5,'None',1};

%%

lgd = cell(nmethod*length(sigma)+2,1);
h1 = zeros(nmethod*length(sigma)+2,1);
h2 = zeros(nmethod*length(sigma)+2,1);
h3 = zeros(nmethod*length(sigma)+2,1);
h4 = zeros(nmethod*length(sigma)+2,1);

fig = figure('Units','inches', 'Position',[0 0 7.5 5.4],'PaperPositionMode','auto');
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',9, 'FontName','Times');

load([savepath, 'objective_results_sigma=0',str_add,'.mat']);

subplot(221);
h1(1) = semilogx(rt60, mean(ops_ideal,2));grid on;hold on;
subplot(222);
h2(1) = semilogx(rt60, mean(tps_ideal,2));grid on;hold on;
subplot(223);
h3(1) = semilogx(rt60, mean(ips_ideal,2));grid on;hold on;
subplot(224);
h4(1) = semilogx(rt60, mean(aps_ideal,2));grid on;hold on;
lgd{1} = 'Ideal TF';

subplot(221);
h1(2) = semilogx(rt60, mean(ops_mwf,2));grid on;hold on;
subplot(222);
h2(2) = semilogx(rt60, mean(tps_mwf,2));grid on;hold on;
subplot(223);
h3(2) = semilogx(rt60, mean(ips_mwf,2));grid on;hold on;
subplot(224);
h4(2) = semilogx(rt60, mean(aps_mwf,2));grid on;hold on;
lgd{2} = 'MCWF';
str_add = ['_', calib_type, '_', method{1}]; 




for i = 1:nmethod
    
    str_add = ['_', calib_type, '_', method{i}];
    
    for k = 1:nsigma

        resname = ['objective_results_sigma=',num2str(sigma(k)),str_add,'.mat'];
        load([savepath, resname]);      
        lgd{(i-1)*nsigma+(k+2)} = [upper(method{i}),' $\sigma=$',num2str(sigma(k))];

        fig;
        newsubplot(221,'','OPS');
        h1((i-1)*nsigma+(k+2)) = semilogx(rt60, mean(ops_mle,2));grid on;hold on;
        ylim([0 200]);xticks(rt60);yticks([0:50:200]);xlim([0, 10]);
    %         xlabel('Source-mic distance (m)');
    %         title('Overall Perceptual Score');

        newsubplot(222,'','TPS');
        h2((i-1)*nsigma+(k+2)) = semilogx(rt60, mean(tps_mle,2));grid on;hold on;
    %         title('Target-related Perceptual Score');
        ylim([0 100]);xticks(rt60);xlim([0, 10]);

        newsubplot(223,'','IPS');
        h3((i-1)*nsigma+(k+2)) = semilogx(rt60, mean(ips_mle,2));grid on;hold on;
    %         title('interference-related Perceptual Score');
        ylim([0 100]);xticks(rt60);xlim([0, 10]);

        newsubplot(224,'','APS');
        h4((i-1)*nsigma+(k+2)) = semilogx(rt60, mean(aps_mle,2));grid on;hold on;
    %         title('Artifact-related Perceptual Score');
        ylim([0 100]);;xticks(rt60);xlim([0, 10]);


        % Give common xlabel, ylabel and title to your figure
        han=axes(fig,'visible','off'); 
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='off';
        xlabel(han,'RT60 (s)','interpreter','latex','FontSize',9);

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

   
% print(['../../figures/',calib_type,'/var_rt60_sub_scrores','.eps'], '-depsc');
