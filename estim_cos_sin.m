% Least-squares estimation of the amplitudes of a cosine and sine 
% model with given frequency
% Input parameters:
% signal: the vector containing the data points
% t: the vector containing the corresponding time samples
% f0 : the frequency of the sine wave
% Output parameters:
% alpha_cos : the amplitude of the cosine part
% beta_sin : the amplitude of the sine part
% in a model alpha_cos*cos(2*pi*f0) + beta_sin*sin(2*pi*f0).



function [alpha_cos,beta_sin] = estim_cos_sin(signal,t,f0)

ind_NZ = find(signal~=0);
t_NZ = t(ind_NZ); t_NZ = t_NZ(:);
signal = signal(ind_NZ); signal = signal(:);

R = [cos(2*pi*f0*t_NZ), sin(2*pi*f0*t_NZ)];
amp_cossin = (R'*R)\(R'*signal);

alpha_cos = amp_cossin(1);
beta_sin = amp_cossin(2);
