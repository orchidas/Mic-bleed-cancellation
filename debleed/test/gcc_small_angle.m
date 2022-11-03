%GCC-PHAT for FIR filters with small angle approximationg

close all;
lags = (-10:0.1:10).';
f = 0.01*(1i.^lags)./lags;

figure;
subplot(211);
plot(lags,[real(f), imag(f)]);ylim([-1,1]);
xlabel('Lags');ylabel('GCC-PHAT for small angles');  grid on;
subplot(212);
plot(lags, abs(f));grid on;