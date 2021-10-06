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
%     mean_y = repmat(mean(y,1), [N,1]);
    Rxx = (y*y')/T;
    % interference correlation matrix, NxN
%     mean_int = repmat(mean(interf,1),[N,1]);
    Rnn = (interf*interf')/T;
  

%     % invert interference matrix
%     [U, S, Ut] = svd(Rnn);
%     Rnn_inv = Ut(:,1:M)*diag(1./diag(S(1:M,1:M)))*U(:,1:M)';

    % GEVD
    [X, L] = eig(Rxx, Rnn);  
    %sort eigenvalues and vectors
    [lambda, order] = sort(diag(L),'descend'); 
    X = X(:, order); 
    
%     %scale eigenvectors, such that X'*Rnn*X = I
%     sc = sqrt(diag(X'*Rnn*X));
%     X = X*diag(1./sc);
%     %new eigenvalues according to scaled eigenvector
%     lambda_norm = diag(X'*Rxx*X);
        
    %inverse filter
    E = complex([eye(M);zeros(N-M,M)]);
    D = [(lambda(1:M)-1)./(lambda(1:M) + mu - 1); zeros(N-M,1)];
    W = (X * diag(D) / X) * E;    
%     L_sort = diag(lambda(1:M));
%     W = (X / (L_sort+(mu-1)*eye(M)) * (L_sort-eye(M)) / X)*E;
    src = W'*y; % separated source



end

