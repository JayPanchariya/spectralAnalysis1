% Computes the Fourier transform and shifts the negative frequencies so
% that the frequency axis is between -Fs/2 and Fs/2

% x input signal
% Nf number of frequencies
% Fs sampling frequency (for defining the frequency axis)

function [freq_centered,FT_centered] = compute_fft_and_shift(x,Nf,Fs)

FTx = fft(x,Nf);
freq = 0:1/Nf*Fs:(Nf-1)/Nf*Fs;

freq_centered = freq - Fs/2;
FT_centered = fftshift(FTx);