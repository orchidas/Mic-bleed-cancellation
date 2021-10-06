function [s_hat, H_opt, H0] = debleed(xmic, Nsrc, fs, frameSize, hopSize, fftSize, method, sigma, calib_flag, varargin)

%% 
% Function to reduce microphone bleed and estimate MIMO transfer function'
%
% Outputs -
% s_hat        Matrix of separated sources (one source in each column)
% H_opt        NxM Optimized transfer function matrix
%
% Inputs -
% xmic         Matrix of mic signals. Each mic signal is a column
% Nsrc         Number of sources we want to isolate
% frameSize    STFT frameSize
% hopSize      STFT hop size
% method       Maximum likelihood ('mle'), MAP ('map'), MWF with GEVD
% sigma        Regularization hyperparameter
% calib_flag   How to estimate initial transfer function (0 or 1)
% calib_mat    Calibration files (if available)
%
% Author       Orchisama Das, 2020
%%

    if calib_flag
        calib_mat = varargin{1};
        calib_type = varargin{2};
%         Sim = varargin{3};
        
    end
    
    
    Nmic = size(xmic,2);   %number of mics
    win = hann(frameSize);
    X = [];      %frequency-domain observation matrix of size Nframes x nbins x Nmics

    % get observation matrix
    for i = 1:Nmic
        M = get_stft_from_audio(xmic(:,i), frameSize, hopSize, fftSize, win);
        if i == 1
            X = zeros(size(M,1),size(M,2),Nmic);
            X(:,:,1) = M;
        else
            X(:,:,i) = M;
        end
    end
    
    Nframes = size(X,1);
    nbins = size(X,2);


    %% find initial estimate of transfer function from each source to each mic
    
    H0 = zeros(nbins,Nmic,Nsrc);  %initial estimate of transfer function matrix from each source to each mic (nbins x nsources x nmics)
    freqaxis = linspace(-fs/2,fs/2,nbins);

    if ~calib_flag
        % if calibration has not been done,
        % estimate solo excerpts
        % a time frame contains a solo instrument if the mic closest to that
        % instrument is significantly louder than all other mics.
        % then, average tf estimates over all time-frames of solo excerpts

        %% find short time energy 

        energy = zeros(Nframes,Nmic); %short time energy in dB

        for mic = 1:Nmic
            for frameno = 1:Nframes
                energy(frameno,mic) = mag2db(sum(abs(X(frameno,:,mic).^2))/nbins);    %parseval's theorem
            end
        end


        %% find solo excerpts

        uthres = 20;    %these numbers need to be tinkered with
        lthres = -20;
        
        %loop over all instruments 
        for mic = 1:Nmic
            nsolo_flag = 0;
            solo_flag = false;
            othermics = setdiff(1:Nmic,mic);

            %loop over all time frames
            for frameno = 1:Nframes    
                if energy(frameno,mic) > uthres && energy(frameno, othermics) < lthres
                    solo_flag = true;
                    nsolo_flag = nsolo_flag+1;
                    H0(:,:,mic) = H0(:,:,mic) + estimate_initial_TFs(reshape(X(frameno,:,:), nbins,Nmic), mic); % need to change this line
                end  
            end

            %if there are any solo excerpts
            if solo_flag
                H0(:,:,mic) = H0(:,:,mic)/nsolo_flag;
                figure(mic);
                semilogx(freqaxis, fftshift(mag2db(abs(H0(:,:,mic)).^2),1));grid on;
                xlabel('Frequency (hz)'); ylabel('Magnitude (dB)');
                xlim([20,20000]);
            else
                warning(['No solo excerpts found in instrument ', num2str(mic), '!']);
            end
        end



    else
        %% if calibration has been done, then the initial estimate is handed to us
        
        H0 = zeros(nbins,Nmic,Nsrc);  %initial estimate of transfer function matrix from each source to each mic (nbins x nsources x nmics)
        legendcell = strings(Nsrc,Nmic);
        frac = Nmic/Nsrc;
        chan_L = 2^7;
        nreflect = 20;
        
        %loop over sources
        for nsrc = 1:Nsrc   
            
            closest_mic = (nsrc-1)*frac+1;
            den = calib_mat.src{nsrc}.mic(:,closest_mic); %pick the closest mic (irrelevant for my setup)
            DEN = get_stft_from_audio(den, frameSize, hopSize, fftSize, win);   
            H0(:,closest_mic,nsrc) = ones(nbins,1); %all the close mics have to have a flat response
            legendcell(nsrc, closest_mic) = ['$H_{', num2str(nsrc), num2str(closest_mic),'}$'];
            
            interf_data_length =  min(size(calib_mat.src{nsrc}.mic,1));   %length of data for calculating interference correlation matrix
            interf = zeros(interf_data_length,Nmic);
            
            %loop over mics
             for nmic = 1:Nmic
                    if  nmic ~= closest_mic

                        num = calib_mat.src{nsrc}.mic(:,nmic);
                        if strcmp(calib_type,'spec-ratio')
                             
                             NUM = get_stft_from_audio(num, frameSize, hopSize, fftSize, win);
%                              num_energy_dB = get_signal_energy(NUM);  %only frames that have significant leakage energy should be used for calibration
%                              frames_to_consider = find(num_energy_dB >= max(num_energy_dB)-3);   %half the maximum loudness
                             frames_to_consider = 1:size(NUM,1); % for drums, consider all frames
                             H0(:,nmic,nsrc) = mean(NUM(frames_to_consider,:)./DEN(frames_to_consider,:),1);
                        
                        elseif strcmp(calib_type, 'gcc-phat')
                            H0(:,nmic,nsrc) = estimate_tf_gain_delay(den, num, fs, nreflect, nbins, chan_L);
                        end
                        
                         interf(:,nmic) = interf(:, nmic) + calib_mat.src{nsrc}.mic(1:interf_data_length,nmic);
                         legendcell(nsrc, nmic) = ['$H_{', num2str(nsrc), num2str(nmic),'}$'];

                    else
                        continue;
                    end                                                                              
             end
            
             if strcmp(calib_type, 'bci')
                % estimate transfer function matrix with BCI
                 [~, H0(:,:,nsrc)] = estimate_tf_bci(calib_mat.src{nsrc}.mic,fs,Nmic,chan_L,nbins);
             end
             
%             figure(nsrc);clf;             
%             H0_2D = reshape(H0(:,:,nsrc),[nbins,Nmic]);
%             h0 = ifft(H0_2D,nbins,1);
%             hideal_2D = reshape(Sim.h_ideal(:,:,nsrc), [chan_L, Nmic]);
%             Hideal_2D = fft(hideal_2D, nbins);
%             
%             for nmic = 1:Nmic
%                 subplot(Nmic,1,nmic);
% %                 stem([[hideal_2D(:,nmic); zeros(nbins-chan_L,1)],h0(:,nmic)]);grid on; hold on;
% %                 xlim([0, chan_L]);ylim([-1 1]);
%                 semilogx(freqaxis,mag2db(abs([Hideal_2D(:,nmic),H0_2D(:,nmic)])));grid on; hold on;
%                 xlim([0, 20000]);               
%             end
%             sgtitle(['Source ', num2str(nsrc)]); 
%             xlabel('Samples'); ylabel('Amplitude');
%             legend('Actual','Spectral Ratio');
            
            
%             semilogx(freqaxis, fftshift(mag2db(abs(H0_2D).^2),1));grid on;hold on;
%             xlabel('Frequency (hz)'); ylabel('Magnitude (dB)');title('Initial transfer function');
%             xlim([20,20000]);
%             legend(legendcell(nsrc,:), 'interpreter','latex');drawnow;

        end


    end
    
    
    %% batch processing over all frames
    
    %% Maximum likelihood estimate
    
    if strcmp(method, 'mle')
        
        H_opt = H0;
        S = zeros(Nframes, nbins, Nsrc);

        % loop over frequencies
        parfor k = 1:nbins/2

            x = reshape(X(:,k,:), [Nframes,Nmic]).';
            Hopt_2D = reshape(H_opt(k,:,:),[Nmic,Nsrc]);

            [S_hat, H_hat] = mle_closed_form_batch(Nsrc, x, Hopt_2D, sigma);
            H_opt(k,:,:) = H_hat;
            S(:,k,:) = reshape(S_hat.', [Nframes,1,Nsrc]);
            disp(['Frequency bin ' ,num2str(k),' has been processed']);

        end
        temp = S(:,1:nbins/2,:);
        S(:,nbins/2+1:end,:) = fliplr(conj(temp)); %conjugate symmetric spectrum 

  
    %% GEVD based Multichannel Wiener filter
    
%     elseif (strcmp(method, 'gevd_mwf'))
%         
%         S = zeros(Nframes, nbins, Nsrc);
%         mu = sigma;   %more weight given to speech distortion
%         for nmic = 1:Nmic
%             int = get_stft_from_audio(interf(:,nmic), frameSize, hopSize, fftSize, win);
%             if nmic == 1
%                 Interf = zeros(size(int,1),size(int,2),Nmic);
%                 Interf(:,:,1) = int;
%             else
%                 Interf(:,:,nmic) = int;
%             end
%         end
% 
%         parfor k = 1:nbins/2
%             x = reshape(X(:,k,:), [Nframes,Nmic]).';
%             interf_sig = reshape(Interf(:,k,:), [size(int,1),Nmic]).';
%             S_hat = gevd_mwf(x,interf_sig,Nmic,Nsrc,mu);
%             S(:,k,:) = reshape(S_hat.', [Nframes,1,Nsrc]);
%             disp(['Frequency bin ' ,num2str(k),' has been processed']);        
%         end
%         
%         S(:,nbins/2+1:end,:) = S(:,1:nbins/2,:); %symmetric spectrum 

    
    %% Maximum aposteriori probability
    
    elseif (strcmp(method,'map'))
        
        H_opt = H0;
        sigma_w = 0.1;
        S = zeros(Nframes, nbins, Nsrc);
        
%        sigma_nu = 0.1;  %variance in \tilde{H} measurement
%        [mu_s, P_s, P_H] = estimate_apriori_statistics(X, H0, sigma_w, sigma_nu, P_Htilde);
              
        for nmic = 1:Nmic
            int = get_stft_from_audio(interf(:,nmic), frameSize, hopSize, fftSize, win);
            if nmic == 1
                Interf = zeros(size(int,1),size(int,2),Nmic);
                Interf(:,:,1) = int;
            else
                Interf(:,:,nmic) = int;
            end
        end
        
       % loop over frequencies
        parfor k = 1:nbins/2
    
            x = reshape(X(:,k,:), [Nframes,Nmic]).';
            Hcur_2D = reshape(H_opt(k,:,:),[Nmic,Nsrc]);
            interf_sig = reshape(Interf(:,k,:), [size(int,1),Nmic]).';
            [Rss_inv, mu_s] = estimate_source_statistics(x, interf_sig, Nsrc);

%             Ps_2D = reshape(P_s(k,:,:), [Nsrc,Nsrc]);
%             PH_2D = reshape(P_H(k,:,:), [Nsrc,Nmic]);
%             [S_hat, H_hat] = map_closed_form_batch(x, Hcur_2D, sigma_w, mu_s(:,k), Ps_2D, PH_2D);
            
            [S_hat, H_hat] = map_closed_form_batch(x,Hcur_2D,sigma,Rss_inv,mu_s,Nsrc,sigma_w);
            H_opt(k,:,:) = H_hat;
            S(:,k,:) = reshape(S_hat.', [Nframes,1,Nsrc]);
            disp(['Frequency bin ' ,num2str(k),' has been processed']);
    
        end
        temp = S(:,1:nbins/2,:);
        S(:,nbins/2+1:end,:) = fliplr(conj(temp)); %conjugate symmetric spectrum 
        
        
    end
    

%% convert STFT to time domain 
s_hat = [];
for src = 1:Nsrc
    s_hat = [s_hat, get_audio_from_stft(S(:,:,src), hopSize)];
end

       
end
    
    
        


