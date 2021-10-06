function [strikeStart, strikeEnd, smooth_level] = onset_detector(snd,fs,tau_a,tau_f,frac_tc,fname)
%% 
% Onset detection with leaky integrator level detector. Takes a drum calibration audio
% as input and returns onset and offset times of the drum hits
% tau_a - attack time (in s)
% tau_f - decay time (in s)
% frac - fraction of time constant we are looking for
%%

N = length(snd);
t = 0:1/fs:(N-1)/fs;
[level,smooth_level] = leaky_level_detector(snd,fs,tau_a,tau_f); %get signal envelope

[pks, onsetSamp] = findpeaks(smooth_level,'MinPeakWidth',0.1*fs, ...
    'MinPeakHeight',0.5*mean(smooth_level));  %find peaks from smoothed level
numOnsets = length(onsetSamp);

%% find offsets

offsetSamp = zeros(numOnsets-1,1);
for k = 1:numOnsets-1
    level_between_onsets = smooth_level(onsetSamp(k)+1:onsetSamp(k+1));
    val_tau = level_between_onsets(1)/exp(1/frac_tc);    %value at frac*time constant
    tau = frac_tc*find(level_between_onsets <= val_tau,1);   %estimate time constant
    offsetSamp(k) = onsetSamp(k) + min(3*tau, length(level_between_onsets)...
        -(0.2*fs)); %extract three time constants' worth
end


%% extract loudest strike

[~, maxOnsetInd] = max(pks(1:end-1));
strikeStart = onsetSamp(maxOnsetInd);
strikeEnd = offsetSamp(maxOnsetInd);


%% plot

% figure;
figure('Units','inches', 'Position',[0 0 3.3 2.4],'PaperPositionMode','auto');
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');

plot(t, snd);hold on; 
plot(t, level, '-r', 'Linewidth',0.8);hold on;

xline((strikeStart-1)/fs, '--k');hold on;
xline((strikeEnd - 1)/fs,'--m'); hold off;
legend('Signal','Level','Onset','Offset','FontSize',6);
% for k = 1:numOnsets
%     xline((onsetSamp(k)-1)/fs,'--k');hold on;
%     if k < numOnsets
%         xline((offsetSamp(k)-1)/fs,'--m');hold on;
%     end
% end
% hold off;

grid on;
axis tight;
xlabel('Time(s)');ylabel('Amplitude');
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');
print(['../figures/', fname,'_level.eps'], '-depsc');

end

