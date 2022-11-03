function [peak_pos] = estimate_mode_locations(H_ir, nsrc, win_len)
%% Estimate mode locations from frequency response 
% nsrc- use any non-close mics to estimate peaks, and find common peak
% locations that are within \pm win_len of peak location

nfft = size(H_ir,1);
Nmics = size(H_ir,2);
maxpeaks = nfft/4;
mode_idx = [];

for i = 1:Nmics
    %close mic measurement might mess things up
    if i == nsrc
        continue;
    else    
        [peaks,idx] = findpeaks(abs(H_ir(1:nfft/2,i,nsrc)), 'MinPeakHeight', ...
        max(abs(H_ir(1:nfft/2,1,2))/sqrt(8)), 'NPeaks',maxpeaks);
        % peaks should be common to all measurements
        if isempty(mode_idx)
            mode_idx = idx;
        else
            mode_idx = ismember(mode_idx,idx+(-win_len:win_len));
        end
    end
end

fbins = 1:nfft/2;
peak_pos = fbins(mode_idx);


% figure;
% plot(fbins, 20*log10(abs(H_ir(1:nfft/2,:,nsrc))));hold on;grid on;
% for k = 1:length(peak_pos)
%     xline(peak_pos(k), '-.k'); hold on;
% end
% hold off;
% ylim([-80,0]);
% set(gca,'xscale','log');

end


