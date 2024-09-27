function [hmin] = minimum_phase_IR(Hmag)
%Returns minimum phase impulse response

% Log-transfer function magnitude
ltf = log(Hmag + eps); % Analytic signal
as = hilbert(ltf);
% Conjugated analytic signal
cas = conj(as);
% Compute minimum -phase impulse response
hmin = ifft(exp(cas),'symmetric');
end

