function [h_rir, delays] = closest_viable_rir(g_ir, lambda, beta, t60, fs, min_room_dim, varargin)
%% find closest viable RIR to given IR matrix
% g_ir - nsampx Nextramic x Nsrc matrix of IRs
% lambda - how much to weigh sparsity
% beta - %window length for energy envelope calculation
% (we must be sure that h(t+\tau) < h(t)
%%

p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'plot_data',0);
parse(p,varargin{:});
plot_data = p.Results.plot_data;

fig = figure();
[nsamp, Nmic, Nsrc] = size(g_ir);
h_rir = zeros(nsamp, Nmic, Nsrc);
delays = zeros(Nmic, Nsrc);
c = 343;
%direct path cannot start or end after this
min_peak_time = round(5e-4*fs); % corresponding to a distance of 0.17m from any source
max_peak_time = round(min_room_dim/c * fs);  %mic cannot go beyond room walls
%time constant for correct decay envelope
time_const = t60/log(1000);

for i = 1:Nmic
    for j = 1:Nsrc
        g = reshape(g_ir(:,i,j), [nsamp,1]);
        %optimize for magnitude only
        h = cvx_optimize_rir(g,lambda, beta, time_const, fs);
        %add sign
        h = h .*sign(g);
        h_rir(:,i,j) = h;
        [max_peak_val, max_peak_loc] = max(abs(h(min_peak_time:max_peak_time)));
%         max_dev = max(diff(abs(h(min_peak_time:max_peak_time))));
%         [~, possible_dp] = findpeaks(abs(h(min_peak_time:max_peak_time)),...
%             'MinPeakHeight', max_peak_val/sqrt(2));
        delays(i,j) = max_peak_loc + min_peak_time - 1;
        if plot_data
            fig;clf;
            subplot(211);plot(g(1:end-beta));title('exact solution');xlim([min_peak_time, nsamp/2]);
            subplot(212);plot(h(1:end-beta));hold on;title('viable RIR');xlim([min_peak_time, nsamp/2]);
            xline(delays(i,j), '-.r');hold off;
            xlabel('Samples'); ylabel('Amplitude');
        end
    end
end

end

