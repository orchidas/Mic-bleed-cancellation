function [X, ax, cb, wx] = ftgram(x, fs, typename, varargin)
% FTGRAM - compute, plot short-time Fourier rransform
%
% X = ftgram(x, FS, TYPE, VARARGIN) returns X, the short-time Fourier
% transform of the input x, computed using stft().  The STFT is plotted using
% the sampling rate FS in Hz according to the string TYPE and parameters
% specified in name, value pairs.  The variable TYPE can be 'rir', 'music' or
% 'speech'; it sets the defaults for the spectrogram and plotting parameters
% described below.
%
% NAME = [RIR_DEFAULT MUSIC_DEFAULT SPEECH_DEFAULT];   %% DESCRIPTION, UNITS
%
% STFT parameters
% 'nbins' = [512 2048 512];  %% dft half length, bins
% 'nskip' = nbins/2;  %% stft hop size, samples
% 
% spectrogram image axes
% 'dbrange' = [80 60 60]; %% gram dynamic range, dB
% 'normalize' = [true true true]; %% normalize spectrogram, indicator
% 'logf' = [true true false]; %% logarithmic frequency axis, indicator
% 'logt' = [true false false];    %% logarithmic time axis, indicator
% 'ms' = [true false false];   %% time/frequency axis ms/kHz scaling, indicator
% 
% waveform onset trimming
% 'trim' = [true false false]; %% trim waveform onset, indicator
% 'preroll' = 10; %% onset zeropad duration, milliseconds
% 'onsetlevel' = 1e-2;	%% onset level, fraction
% 
% waveform plot parameters
% 'waveform' = [true false false];    %% plot waveform, indicator
% 'tanhflag' = [true false false];    %% hyperbolic tangent saturation, indicator
% 'tanhbeta' = 5; %% hyperbolic tangent saturation parameter, ratio
%
% [~, AX, CB, WX] = ftgram() returns AX, the plot axes, one for each spectrogram
% channel, the first being the waveform plot, if present, and CB, the associated
% color bars, and WX, the waveform lines handle.
%
% See Also: STFT, IRGRAM.

% Created:  2-Feb-2011.
% Revised:  2-Feb-2011, MJW w/JSA, v1.
% Revised: 14-Mar-2012, JSA, v2 - interface, windowing reworked.
% Revised: 11-Jun-2012, JSA, v3 - nbins default changed for 'music' setting.
% Revised: 22-Aug-2013, JSA, v4 - seismic frequency range flag added.
% Revised:  8-Feb-2015, JSA, v5 - colorbar handles returned.
% Revised: 24-May-2017, JSA, v6 - Matlab 2016b compatibility added; XTick fixed.
% Revised: 24-May-2017, JSA, v7 - 'normalize' spectrogram normalization added.
% Revised:  1-Dec-2017, EKCD,v8 - colorbar label position fixed.
% Revised: 21-Dec-2017, JSA, v9 - surf image cropping "bug" accomodated.
% Version: v9.


%% initialization, parse input
% initialize type defaults
nbins_default = [512 2048 512];    %% dft half length, bins
% nskip_default = nbins/2;      %% stft hop size, samples

dbrange_default = [80 60 60];   %% gram dynamic range, dB
normalize_default = [true true true];   %% normalize spectrogram, indicator
logf_default = [true true false];   %% logarithmic frequency axis, indicator
logt_default = [true false false];  %% logarithmic time axis, indicator
logtmin_default = [10 10 10];  %% logarithmic time axis offset, milliseconds
ms_default = [true true false];    %% time/frequency axis ms/kHz scaling, indicator

seismic_default = false;    %% seismic frequency axis, indicator

trim_default = [true false false];  %% trim waveform onset, indicator
preroll_default = 10*[1 1 1];   %% onset zeropad duration, milliseconds
onsetlevel_default = 1e-2*[1 1 1];  %% onset level, fraction

waveform_default = [true false false];  %% plot waveform, indicator
tanhflag_default = [true false false];  %% hyperbolic tangent saturation, indicator
tanhbeta_default = 2*[1 1 1];   %% hyperbolic tangent saturation parameter, ratio

% set signal type
switch typename,
    case 'rir',
        % room impulse response
        type = 1;
    case 'music',
        % music input
        type = 2;
    case 'speech',
        % speech input
        type = 3;
    otherwise,
        % music default
        type = 2;
end;

% parse input
p = inputParser;

% waveform, sampling rate
p.addRequired('x', @(x)isnumeric(x));   %% waveform, signal matrix
p.addRequired('fs', @(x)isnumeric(x) && x>0);   %% sampling rate, Hz
p.addRequired('typename', @(x)ischar(x));   %% sampling rate, Hz

% stft parameters
p.addParameter('nbins', nbins_default(type), @(x)isnumeric(x));
p.addParameter('nskip', 0, @(x)isnumeric(x));

% spectrogram image axes
p.addParameter('dbrange', dbrange_default(type), @(x)isnumeric(x));
p.addParameter('normalize', normalize_default(type), @(x)islogical(x));
p.addParameter('logf', logf_default(type), @(x)islogical(x));
p.addParameter('logt', logt_default(type), @(x)islogical(x));
p.addParameter('logtmin', logtmin_default(type), @(x)isnumeric(x));
p.addParameter('ms', ms_default(type), @(x)islogical(x));

p.addParameter('seismic', seismic_default, @(x)islogical(x));


% waveform onset trimming
p.addParameter('trim', trim_default(type), @(x)islogical(x));
p.addParameter('preroll', preroll_default(type), @(x)isnumeric(x));
p.addParameter('onsetlevel', onsetlevel_default(type), @(x)isnumeric(x));

% waveform plot parameters
p.addParameter('waveform', waveform_default(type), @(x)islogical(x));
p.addParameter('tanhflag', tanhflag_default(type), @(x)islogical(x));
p.addParameter('tanhbeta', tanhbeta_default(type), @(x)isnumeric(x));

p.parse(x, fs, typename, varargin{:});

% assign variables
nbins = p.Results.nbins;
if (p.Results.nskip > 0);
    nskip = p.Results.nskip;
else,
    nskip = nbins/2;
end;

dbrange = p.Results.dbrange;
normalize = p.Results.normalize;
logf = p.Results.logf;
logt = p.Results.logt;
logtmin = p.Results.logtmin;
ms = p.Results.ms;

seismic = p.Results.seismic;

trim = p.Results.trim;
preroll = p.Results.preroll;
onsetlevel = p.Results.onsetlevel;

waveform = p.Results.waveform;
tanhflag = p.Results.tanhflag;
beta = p.Results.tanhbeta;

% plot font sizes
axislabelsize = 8;
ticklabelsize = 8;
cbarlabelsize = 7;
% beta = 3;


%% condition input, form spectrogram

% find input signal size
nsamp = size(x,1);
if (nsamp == 1),
    % make x a column
    x = x(:);
    [nsamp, nc] = size(x);
end;
[nsamp, nc] = size(x);  %% signal length, samples; channel count, channels

% trim waveform onset
if trim,
    istart = find(mean(abs(x)/max(max(abs(x))),2) > onsetlevel, 1) - round(preroll*fs/1000);
    x = x(max(1,istart):end,:);
end;
nsamp = size(x,1);

% compute, normalize spectrogram
X = stft(x, nbins, nskip);
nframes = size(X,2)/nc;

if normalize,
    y = x/max(max(abs(x)));
    Y = 20*log10(abs(X)/max(max(abs(X)))+eps);
else,
    y = x/max(max(abs(x)));
    Y = 20*log10(abs(X));
end;


%% plot waveform

if waveform,
    % define axis
    ax = subplot(max(nc+1,3),1,1);

    % tanh scaling
    if tanhflag,
        y = tanh(beta*y)/tanh(beta);
    end;

    % define time axis
%     y_off = [zeros(round(0.01*fs),1);y];    %Orchi - offset RIR by 10ms
    nsamp = length(y);
    t = [0:nsamp-1]/fs;
    tscale = (~ms)*1 + ms*1000;
    

    % plot waveform, label time axes
    if logt,
%         wx = semilogx(tscale*(t+logtmin/1000), y_off); grid;

        wx = semilogx(tscale*(t + logtmin/1000), y); grid;
        
        xrange = tscale*[logtmin/1000 logtmin/1000+nsamp/fs];
        xlim(xrange);

        temp = kron(10.^(-3:6)', [1 2 5]');
        index = find((temp >= xrange(1)) & (temp <= xrange(2)));
        nxtk = length(index);
        xtk = temp(index);
        set(gca, 'XTick', xtk);

        xtklabel = cell(length(xtk),1);
        for tk = (1:nxtk),
            xtklabel{tk} = sprintf('%g', xtk(tk));
        end;
        set(gca,'XTickLabel',xtklabel);

    else,
        plot(tscale*t, y); grid;
        xlim(tscale*[0 nsamp/fs]);
    end;

    % y-axis labels
    if tanhflag,
        ylabel('amplitude, dB', 'FontSize', axislabelsize);
%          set(gca,'YTick', sort(kron([-1 1], tanh(10.^-([0:10:30]/20)*beta))));
%          set(gca,'YTickLabel',['  0'; '-10'; '-20'; '-30'; '-30'; '-20'; '-10'; '  0'])
%          set(gca,'YTick', sort(kron([-1 1], tanh(10.^-([0:15:30]/20)*beta))));
%          set(gca,'YTickLabel',['  0'; '-15';  '-30'; '-30';  '-15'; '  0'])
         set(gca,'YTick', sort(kron([-1 1], tanh(10.^-([0:20:20]/20)*beta))));
         set(gca,'YTickLabel',['  0'; '-20';  '-20'; '  0'])
        ylim(tanh(10^(1/20)*beta)*[-1 1]);
    else,
        ylim([-1 1]);
        ylabel('amplitude', 'FontSize', axislabelsize);
    end;

end;


%% plot spectrogram

% define time, frequency, energy axes
tscale = (~ms)*1 + ms*1000;

np = ceil(nbins/(2*nskip));
nq = ceil((nsamp-(nframes-1)*nskip-nbins/2)/nskip);
t = tscale * ([-np ((0.5-np):nframes-1+nq-0.5) nframes-1+nq]*nskip + nbins/2)/fs;
f = 1/tscale * fs/2*[0 (0.5:nbins-0.5) nbins]/nbins;
f(1) = eps;

cb = zeros(1,nc);

% loop through specgtrograms
for s = [1:nc],
     set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');

    
    % get axes
    if (nc == 1),
        if waveform,
            ax(2) = subplot(3,1,[2 3]);
        else,
            ax(1) = subplot(1,1,1);
        end;
    else,
        ax(waveform+s) = subplot(nc+waveform,1,s+waveform);
    end;

    % display spectrogram
    offset = logt * tscale*logtmin/1000;
    surf(t+offset, f, Y([(1:nbins+1) nbins+1], (s-1)*nframes + [ones(1,np) [1:nframes] nframes*ones(1,nq+1)]), 'edgecolor', 'none');
    axis tight;
    view(0,90);

    % time, frequency scaling
    if ms,
        if (s == nc),
            xlabel({'time, ms',''}, 'FontSize', axislabelsize);
        end;
        ylabel('frequency, kHz', 'FontSize', axislabelsize);
    else,
        if (s == nc),
            xlabel('time, s', 'FontSize', axislabelsize);
        end;
        ylabel('frequency, Hz', 'FontSize', axislabelsize);
    end;

    % scale, label frequency axis
    if logf,
        set(gca, 'yScale', 'log');
    end;

    % scale, label time ax1s
    if logt,
        set(gca, 'xScale', 'log');

        xrange = tscale*[logtmin/1000 logtmin/1000+nsamp/fs];
        xlim(xrange);

        temp = kron(10.^(-3:6)', [1 2 5]');
        index = find((temp >= xrange(1)) & (temp <= xrange(2)));
        nxtk = length(index);
        xtk = temp(index);
        set(gca, 'XTick', xtk);

        xtklabel = cell(length(xtk),1);
        for tk = (1:nxtk),
    
            xtklabel{tk} = sprintf('%g', xtk(tk));
   
        end;
        set(gca,'XTickLabel',xtklabel);

    end;

    % scale, label frequency axis
    if logf,
        if seismic,
            % seismic frequency axis
            divs = [10 20 50 100 200 500 1000]/tscale;
            set(gca, 'ytickmode', 'manual');
            set(gca, 'ytick', divs);

            ylim([min(divs) min(fs/2,max(divs))]);

        else,
            % audio frequency axis
%             divs = [20 50 100 200 500 1000 2000 5000 10000 20000]/tscale;
            divs = [50 100 500 1000 5000 10000]/tscale;
            set(gca, 'ytickmode', 'manual');
            set(gca, 'ytick', divs);

            ylim([min(divs) max(divs)]);

        end;

    else,
        if seismic,
            % seismic frequency axis
            ylim([0 min(fs/2,1000)]/tscale);

        else,
            % audio frequency axis
%             ylim([0 20000]/tscale);
            divs = [500 2000 5000 8000]/tscale;
            set(gca, 'ytickmode', 'manual');
            set(gca, 'ytick', divs);
            set(gca, 'yticklabels', divs/1000);
            ylim([0 8000]/tscale);

        end;
    end;

    % display color bar
    caxis([-dbrange 0])
    cb(s) = colorbar();
    if waveform*(nc == 1),
        set(cb(s), 'Position', [0.916 0.11 0.015 0.517]);
    elseif (nc == 1),
        set(cb(s), 'Position', [0.916 0.11 0.015 0.815]);
    else,
        temp = get(cb(s), 'Position');
        set(cb(s), 'Position', [0.916 temp(2)+0.004 0.015 temp(4)]);
    end;
    colormap(jet);
%     ylabel(cb(s),'energy, dB', 'FontSize', axislabelsize);
%     set(get(cb(s), 'Label'), 'Position', [3.75 -35.0 0]);

end;

% % adjust axis font sizes
for iax = (1:length(ax)),
    set(ax(iax), 'FontSize', ticklabelsize);
end;

for icb = (1:length(cb)),
    set(cb(icb), 'FontSize', cbarlabelsize);
end;

%link x-axes of plots:
if (waveform + nc-1),
    linkaxes(ax,'x');
end;
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');


end
