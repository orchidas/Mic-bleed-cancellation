function [x] = get_audio_from_stft(X, hopSize)

%convert STFT back to audio signal

    [nFrames, fftSize] = size(X);
    x = zeros(nFrames*hopSize + (fftSize - hopSize), 1);
%     start = 1;

    for i = 1:nFrames
        start = (i-1)*hopSize + 1;
        x_frame = (ifft(X(i,:), fftSize));
        x(start:start + fftSize - 1) = x(start:start + fftSize - 1) + x_frame.';
%         start = start + hopSize;
    end
    x = real(x);

end


