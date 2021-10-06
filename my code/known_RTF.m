function [s_ideal] = known_RTF(xmic, Nsrc, h_ideal, frameSize, hopSize, fftSize)

%% Separate sources when the ideal transfer function from each source to mic is known
% addpath(genpath('../BSIE_toolbox/.'));   %BSI toolbox

Nmic = size(xmic,2);
win = hann(frameSize);
X = [];      %frequency-domain observation matrix of size Nframes x nbins x Nmics

% get observation matrix
for i = 1:Nmic
    M = get_stft_from_audio(xmic(:,i), frameSize, hopSize, fftSize, win);
    if i == 1
        X = zeros(size(M,1),size(M,2),Nmic);
        X(:,:,1) = M;
    else
        X(:,:,i) = M;
    end
end

%% use system equalization - does not work for MIMO systems

% Li = length(h_ideal);
% k = 0;
% Nmic = size(xmic,2);
% 
% gh = lsinvfilt(h_ideal, Li, k);
% g = reshape(gh, [Li, Nsrc, Nmic]);


Nframes = size(X,1);
nbins = size(X,2);
H_ideal = fft(h_ideal, fftSize, 1);
S = zeros(Nframes,nbins,Nsrc);

for k = 1:nbins/2
    x = reshape(X(:,k,:), [Nframes,Nmic]).'; 
    Hideal_2D = reshape(H_ideal(k,:,:), [Nmic, Nsrc]);
    s_hat = Hideal_2D\x;
    S(:,k,:) = reshape(s_hat.', [Nframes,1,Nsrc]);
end
temp = S(:,1:nbins/2,:);
S(:,nbins/2+1:end,:) = fliplr(conj(temp)); %conjugate symmetric spectrum 
% plot(abs(S(1,:,1)));

s_ideal = [];
for src = 1:Nsrc
    s_ideal = [s_ideal,get_audio_from_stft(S(:,:,src), hopSize)];
end



end

