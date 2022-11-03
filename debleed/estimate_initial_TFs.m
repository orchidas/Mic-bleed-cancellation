function[H0] = estimate_initial_TFs(Xsolo, srcn)

%Xsolo - observation (mic) signals STFT for a solo excerpt N x nbins
%srcn - source number
%H0 is a vector of size Nxnbins

%For a solo excerpt,the transfer function is calculated as the spectral ratio from
%the closest mic to all the other mics.

[nbins,N] = size(Xsolo);

H0 = zeros(nbins,N);
eps = 10^-12;   %to avoid division by 0

for i = 1:N
    H0(:,i) = Xsolo(:,srcn)./(Xsolo(:,i) + eps);
end

end