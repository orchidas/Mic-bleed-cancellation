function [tau] = gcc_phat(x1, x2, fs, maxpeaks)
%% TDOA estimation with GCC-PHAT 

%%

l1 = length(x1);
l2 = length(x2);
M = l1+l2-1;
nfft = 2^nextpow2(M);

X1 = fft(x1, nfft);
X2 = fft(x2, nfft);
cpsd = conj(X1).*X2;    %cross power spectrum

%cross correlation function
% R = ifft(cpsd./abs(cpsd));  %GCC-PHAT
R = ifft(cpsd./(abs(X1).^2)); %modified GCC

% R = real([R(end-M+2:end); R(1:M)]);
% lags = [(M-1):-1:0, 1:M-1].';
% % remove bias (triangular weighting)
% R = R./(M-lags);

R = real(R(1:M));
lags = (0:M-1)';


% find maximum in GCC function
% the first reflection will always be the largest
[peak, tdoa] = max(abs(R));
[peaks, locs] = findpeaks(abs(R(tdoa-1:end)));
locs = locs + tdoa - 2;

% [peaks, locs] = findpeaks(abs(R));
[peaks, inds] = sort(peaks,'descend');
reflect = locs(inds(1:maxpeaks));


% parabolic interpolation
R_db = mag2db(abs(R));
tau = zeros(maxpeaks,1);
delay_interp = zeros(maxpeaks,1);
peak_interp = zeros(maxpeaks,1);

for i = 1:maxpeaks
    a=R_db(reflect(i)-1);
    b=R_db(reflect(i));
    c=R_db(reflect(i)+1);
    p = 0.5*((a-c)/(a+c-2*b));
    peak_interp(i) = db2mag(b - (0.25*(a-c)*p));
    delay_interp(i) = lags(reflect(i)) + p; %in bins
    tau(i) = delay_interp(i)/fs;
end
    
% figure(1);clf;
% plot(lags/fs, R);grid on;xlabel('Lag (in s)');hold on;
% plot(delay_interp/fs, peak_interp, 'r*');xlim([0, 1]);
% hold off;




end

