function [X] = get_stft_from_audio(x, frameSize, hopSize, fftSize, win)

% A function to get the stft of an audio signal 
% using the Hann window
% x: input signal (row or column vector) - assume real
% fs: sampling rate of x
% frameSize: frame size (in samples) = window length (make it odd)
% hopSize: time between start-times of successive windows (in samples)
% fftSize: sets the zero-padding factor - defaults to frameSize


    N = length(x);
    if(frameSize > N)
        error('Frame size has to be lesser than signal length');
    end
    if(hopSize > frameSize)
        error('Hop size cannot be greater than frame size');
    end

%     make sure window size is odd
%     if(mod(frameSize,2) == 0)
%         frameSize = frameSize+1;
%     end
 
    %convert fftSize to power of 2 
    fftSize = 2^nextpow2(fftSize);

    %find number of frames
    nframes = ceil(N/hopSize);

    %zero pad input signal so that it has nframes exactly
    x = [x; zeros(nframes*hopSize - N + (frameSize - hopSize),1)]; 
    x_frames = zeros(nframes, frameSize);
    X = zeros(nframes, fftSize);
%     start = 1;

    for frame = 1:nframes
        start = (frame-1)*hopSize + 1;
        x_frames(frame, :) = x(start:start+frameSize - 1) .* win;
%         start = start + hopSize;
        X(frame, :) = (fft(x_frames(frame,:), fftSize));
    end

end