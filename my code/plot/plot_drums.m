close all;

readpath = '../../data/Drums_Noah/Saved files/';
load([readpath, 'drums.mat']);

src_name = drums.src_name;
Nmic = drums.Nmic;
Nsrc = drums.Nsrc;
xmic = drums.mic;
fs = 44100;
legendcell = strings(Nsrc,Nmic);

save_file = 0;
sigma = [0 1 100];
method = 'mle';
calib_type = 'bci';
calib_flag = true;
str_add = ['_', calib_type,'_', method]; 
count = 1;



for k = 1:length(sigma)
    
    %% read saved struct
    load([readpath,'separated_drums',str_add,'_sigma=',num2str(sigma(k)),'.mat']);

    
%     H_init = SepDrums.IR_init;
    H_opt = SepDrums.IR_est;
    s_est = SepDrums.mle_source;
    fftSize = size(H_opt,1);
    freqaxis = linspace(-fs/2,fs/2,fftSize);
    
    %H^{-1)x = s
    H_opt_inv = zeros(size(H_opt));   
    for w = 1:fftSize
        H_2D = reshape(H_opt(w,:,:), [Nmic, Nsrc]);
        H_opt_inv(w,:,:) = inv(H_2D);
    end
    
    for nsrc = 4 %1:Nsrc-2

        for nmic = 1:Nmic
            legendcell(nsrc, nmic) = ['$H_{', num2str(nsrc), num2str(nmic),'}$'];         
        end
        
        %% plot initial and optimized transfer functions

        if save_file
            figure('Units','inches', 'Position',[0 0 3.3 2.4],'PaperPositionMode','auto');
            set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',5, 'FontName','Times');
            figure;
            H_2D = reshape(H_opt(:,:,nsrc),[fftSize,Nmic]);
            semilogx(freqaxis, fftshift(mag2db(abs(H_2D).^2),1));grid on;hold on;
            xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
            xlim([20,20000]);
            limsy=get(gca,'YLim');set(gca,'Ylim',[limsy(1) 150]);
            legend(legendcell(nsrc,:), 'interpreter','latex', 'Location',...
                'northwest','NumColumns',2, 'FontSize',6);drawnow;
            title(src_name{nsrc});
            set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',7, 'FontName','Times');
            print(['../figures/',calib_type,'/drums_TF_',method,'_sigma=',...
                num2str(sigma(k)),'.eps'], '-depsc');

        else
        
            figure(2*count-1);clf; 
            for i = 1:2
                subplot(2,1,i);
                if i == 1
                    H_2D = reshape(H_init(:,:,nsrc),[fftSize,Nmic]);
                else
                    H_2D = reshape(H_opt(:,:,nsrc),[fftSize,Nmic]);
                end
                semilogx(freqaxis, fftshift(mag2db(abs(H_2D).^2),1));grid on;hold on;
                xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
                xlim([20,20000]);
            end
            subplot(211);title('Initial transfer functions');
            subplot(212);title('Optimized transfer functions');
            sgtitle([src_name{nsrc}, ', ', method, ', ' ,calib_type, ', $\sigma=$',...
                num2str(sigma(k))],'interpreter','latex');
            legend(legendcell(nsrc,:), 'interpreter','latex');drawnow;
       end

       
       %% plot spectrograms
        [x,~] = audioread([readpath, src_name{nsrc}, '_lnorm.wav']);
        [s,~] = audioread([readpath, upper(method),'/sigma=',num2str(sigma(k)),'/Normalized/',...
            src_name{nsrc},str_add,'_lnorm.wav']);
                                      


        % save spectrograms as individual files
        if save_file

            if sigma(k) == 0
                figure('Units','inches', 'Position',[0 0 2.9 2.5],'PaperPositionMode','auto');
                set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',7, 'FontName','Times');
                ftgram(x,fs, 'music');
                print('../../figures/drums_spec_recorded.eps', '-depsc');
            end
            figure('Units','inches', 'Position',[0 0 2.9 2.5],'PaperPositionMode','auto');
            set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',7, 'FontName','Times');
            ftgram(s,fs,'music');
            print(['../../figures/',calib_type,'/drums_spec_',method,'_sigma=',...
                num2str(sigma(k)),'.eps'], '-depsc');
        else

            figure(2*count);clf;       
            sgtitle(['Recorded (T), Separated (B) ', src_name{nsrc}, ', loudness = -14dB']); 
            ftgram([x s(1:length(x))],fs,'music', 'logf',false);
        end


        if sigma(k) == 0
            s_mat = zeros(length(sigma),length(x));
            s_mat(k,:) = x;
        else
            s_mat(k,:) = s(1:length(x));
        end
        
        count = count+1;
   end
end
      
%% plot recoreded and separated spectrograms in one figure

% figure('Units','inches', 'Position',[0 0 2.9 7.3],'PaperPositionMode','auto');
% set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',7, 'FontName','Times');
% ftgram(s_mat.',fs,'music');
% subplot(311);title('Recorded', 'interpreter','latex');
% subplot(312);title(['Optimized $\sigma=$',num2str(1)],'interpreter','latex');
% subplot(313);title(['Optimized $\sigma=$',num2str(100)],'interpreter','latex');
% print(['../../figures/',calib_type,'/drums_spec_',method,'.eps'], '-depsc');

