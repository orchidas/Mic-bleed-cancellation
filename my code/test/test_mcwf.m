close all, clear all;
addpath('../kokkinis/.');

%%
% INPUTS:
%   m - A matrix containing the microphone signals of size (no.samples x no.microphones)
%
%   WinLen - The length of the window analysis which is also the length of the Wiener filter
%
%   HopSize - The hop size in samples
%
%   Wnd - The window to be used
%
%   beta - The Wiener filter exponent
%
%   Kp - number of averaging frames for Pyy
%
%   wdom - dominant component weight
%
%   wresx - residual component weight for Pxx
%   
%   wresy - residual component weight for Pyy
%
%   psdwemode - if 1 then PSD-WE is applied   
%
% OUTPUTS:
%   out - A matrix containing the processed signals of size (no.samples x
%   no.sources)end
%%

soundpath = '../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Audio/';  %anechoic string quartet recordings
writepath = '../data/26_TU_Berlin_Mozart_Quartett_Konzerthaus/Results/';
fname = {'m1_trim','m2_trim'};
sname = {'Va','Vcl'};
N = length(fname);
xmic =[];


%xmic is input vector - number of samples x number of mics
for i = 1:N
    [m,fs] = audioread([soundpath,fname{i},'.wav']);
    xmic = [xmic m];
end

WinLen = 1024;
HopSize = 512;
Wnd = hann(WinLen);

% from chapter 6, table 6.2 of his thesis
wdom = 0.6;
wresx = 0.25;
wresy = 0.15;

% I am guessing 
Kp = 5;
beta = 1;
psdwemode = 0;

out = mcwf_dmrs(xmic, WinLen, HopSize, Wnd, beta, Kp, wdom, wresx, wresy, psdwemode);

%%
for i = 1:N
    y = out(:,i);
%     soundsc(y,fs);pause;
    audiowrite([writepath,'MCWF_separated_',sname{i},'.wav'],y./max(abs(y)),fs);
end
