function [level,smooth_level] = leaky_level_detector(snd, fs, tau_a, tau_f)
%Level detection with leaky integrator

N = length(snd);
level = zeros(size(snd));

%% level detector

for n = 2:N
    if snd(n) > level(n-1)
        level(n) = level(n-1) + (1-exp(-1/(tau_a*fs)))*(abs(snd(n)) - level(n-1));  
    else
        level(n) = level(n-1) + (1-exp(-1/(tau_f*fs)))*(abs(snd(n)) - level(n-1));   
    end    
end
  
%% fine onsets

smooth_level = smooth(level, 3);    %smooth level with moving average filter
end

