function [hmin] = minimum_phase_IR(Hmag)
%Returns minimum phase impulse response

%method 1
% N = length(Hmag);
% %take the log of the magnitude response
% logmag = log(Hmag);
% argH = -hilbert(logmag);
% hmin = ifft(exp(logmag + 1i*imag(argH)),N);
% hmin = real(hmin);

%method 2
% Log-transfer function magnitude
ltf = log(Hmag + eps); % Analytic signal
as = hilbert(ltf);
% Conjugated analytic signal
cas = conj(as);
% Compute minimum -phase impulse response
hmin = ifft(exp(cas),'symmetric');
end

