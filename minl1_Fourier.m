% Minimisation du critere 0.5||x - Wu||^2 + \mu ||x||_1, où y est un jeu de
% données régulièrement échantillonnées, W est l'opérateur de Fourier
% inverse reconstruisant le spectre x discrétisé sur l'ensemble des fréquences
% axe_freq

% SB, janvier 2014

function [umin,axe_freq] = minl1_Fourier(x,axe_freq,mu_regul)

x = x(:);  %column vector
N = size(x,1); %length of signal
axe_freq = axe_freq(:); % frequency vector
K = size(axe_freq,1); % no of frequency atoms
axe_t = 0:N-1; % Time indices

W = exp(1j*2*pi*axe_t'*axe_freq')/sqrt(N); % W is a matrix where each column is a complex sinusoid of a given frequency in axe_freq.
 
% sparse coefficient initialization
u = zeros(K,1);  % solution
u_old = zeros(K,1); % for convergence check


CV = 0;
n_it = 0;

e = x - W*u; %residual error



while ~CV
    n_it = n_it + 1;

    if ~mod(n_it,50) || (n_it==1)
        fprintf('iteration %g\n',n_it)
        u_old = u;
    end

    % balayage coordonnées
    for k = 1:K
        wk = W(:,k); % k-th column of dictionary
        if u(k)~=0
            e = e + u(k)*wk;
        end
        wkdage = wk'*e;  % Correlation with residual

        % seuillage doux
        if abs(wkdage) < mu_regul
            u(k) = 0;
        else
            temp = (abs(wkdage)-mu_regul)*exp(1i*angle(wkdage));
            u(k) = temp;
            e = e - temp*wk; % Update residual
        end
    end
%Each iteration updates one frequency component based on correlation with the residual (wk'*e).
%If below threshold, set to 0; otherwise, shrink toward zero (soft threshold).
%This is the coordinate descent with soft-thresholding approach.

    % convergence criterian
    CV = norm(u-u_old)/norm(u)< 1e-3 || norm(u)==0;
end


umin = u;

% debiais
ind_NZ = find( u); % indices of non zero coefficient
umin = zeros(size(u));
W_NZ = W(:,ind_NZ);
umin(ind_NZ) = (W_NZ'*W_NZ)\(W_NZ'*x)/sqrt(N)*2; %Refines the solution using least squares only on active (non-zero) frequencies for better amplitude estimation (removes L1 bias).
%plot( abs(umin))
