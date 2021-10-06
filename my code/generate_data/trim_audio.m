function [x_trim] = trim_audio(x,fs,start,stop)

    % start, stop are in samples

    winsize = round(0.005*fs); %window to avoid clicks
    win = hann(2*winsize);
    winl = win(1:winsize);
    winr = win(winsize+1:end);

    x = x(start:stop,:);
    x_trim = [x(1:winsize).*winl; x(winsize+1:end-winsize); ...
            x(end-winsize+1:end).*winr];
    
end

