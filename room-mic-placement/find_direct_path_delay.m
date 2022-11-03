function [delays] = find_direct_path_delay(g_ir, min_room_dim, fs, varargin)
%Simple maximum peak detection in RIR


p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'plot_data',0);
parse(p,varargin{:});
plot_data = p.Results.plot_data;

fig = figure();
[nsamp, Nmic, Nsrc] = size(g_ir);
delays = zeros(Nmic,Nsrc);
c = 343;
%direct path cannot start or end after this
min_peak_time = round(5e-4*fs); % corresponding to a distance of 0.17m from any source
max_peak_time = round(min_room_dim/c * fs);  %mic cannot go beyond room walls


for i = 1:Nmic
    for j = 1:Nsrc
        g = reshape(g_ir(:,i,j), [nsamp,1]);
        [max_peak_val, max_peak_loc] = max(abs(g(min_peak_time:max_peak_time)));
        max_dev = max(diff(abs(g(min_peak_time:max_peak_time))));
        [~, possible_dp] = findpeaks(abs(g(min_peak_time:max_peak_time)),...
            'MinPeakHeight', max_peak_val/sqrt(2), 'Threshold', max_dev*0.5);
        delays(i,j) = possible_dp(1) + min_peak_time - 1;
        if plot_data
            fig;clf;
            subplot(211);plot(g);title('exact solution');xlim([min_peak_time, nsamp/2]);
            xline(delays(i,j), '-.r');hold off;
            xlabel('Samples'); ylabel('Amplitude');
        end
    end
end


        
end

