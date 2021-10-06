function [energy] = get_signal_energy(X)
%X is the STFT of a signal 
%returns energy of each frame of X in dB

    [Nframes,nbins]= size(X);
    energy = zeros(Nframes,1);
    for frameno = 1:Nframes
        energy(frameno) = mag2db(sum(abs(X(frameno,:).^2))/nbins);    %parseval's theorem
    end

end

