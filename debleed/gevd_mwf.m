function [src] = gevd_mwf(y, interf, N, M, mu)
%% 
% GEVD-based MWF for interference cancellation based on
% HASSANI et al.: GEVD-BASED LOW-RANK APPROXIMATION FOR DISTRIBUTED ADAPTIVE 
% NODE-SPECIFIC SIGNAL ESTIMATION
% Inputs:
% y - matrix of microphone signals in STFT domain,NxT - time frames along rows
% interf - interference matrix, MxT - time frames along rows
% N - number of mics
% M - number of sources
% u - distortion vs interference reduction weighting 
%%

    
    if nargin < 5
        mu = 1;
    end
    
    T = size(y,2);
    % mic correlation matrix, NxN
    Rxx = (y*y')/T;
    % interference correlation matrix, NxN
    Rnn = (interf*interf')/T;
  

    % GEVD
    [X, L] = eig(Rxx, Rnn);  
    %sort eigenvalues and vectors
    [lambda, order] = sort(diag(L),'descend'); 
    X = X(:, order); 
    
        
    %inverse filter
    E = complex([eye(M);zeros(N-M,M)]);
    D = [(lambda(1:M)-1)./(lambda(1:M) + mu - 1); zeros(N-M,1)];
    W = (X * diag(D) / X) * E;    
    src = W'*y; % separated source



end

