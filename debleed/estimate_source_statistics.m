function [Rss_inv, mu_s] = estimate_source_statistics(mic, interf, M)
%%
% Estimate source autocorrelation matrix using GEVD
% Inputs:
% mic - mic signal
% interf - interference signal
% M - number of sources
% Output
% Rss_inv - Inverse of source autocorrelation matrix

    %total number of mics and time frames 
    [N,T] = size(mic);
    % mic autocorrelation matrix, NxN 
    Rxx = (mic*mic')/T;
    % interference correlation matrix, NxN
    Rnn = (interf*interf')/T;
    E = [zeros(N-M,M);eye(M)];


    % GEVD
    [X, L] = eig(Rxx, Rnn);  
    % sort eigenvalues and vectors and take top M 
    [lambda, order] = sort(diag(L),'descend');
    % scale eigenvectors so that Rnn = I
    % sc = diag(X'*Rnn*X);    %I think this was found to be identity anyway
                     
    X = X(:, order); 
    phi = 1./(lambda(1:M) - ones(M,1));
    X = X(:,1:M);
    % inverse of source PSD
    Rss_inv = E' * (X * diag(phi) / X) * E; 
    % mean of source signal
    mu_s = mean(mic,2) - mean(interf,2);

end
