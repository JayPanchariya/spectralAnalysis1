close all
load('signals.mat')
% x=[x1,x2,x3,x4];
function[S_MUSIC,freq,ARspectrum,mu]=highResolution(x)
Fs=1;
a=0.1;
mu=a*max(x);
N = length(x);
Nf = 2^nextpow2(N)*16;
freq = (0:Nf-1)/Nf;

P=4;
[a,sigma2] = arcov(x,P);
ARspectrum = sigma2./abs((fft(a,Nf))).^2;
S_MUSIC = pmusic(x,P,freq,Fs);
ARspectrum= ARspectrum(Nf/2+1 : Nf);
ARspectrum(2:end-1) = 2 * ARspectrum(2:end-1); 
S_MUSIC  = S_MUSIC(Nf/2+1 : Nf);
S_MUSIC(2:end-1) = 2 * S_MUSIC(2:end-1);  % Double all but DC and Nyquist
freq   = freq(Nf/2+1 : Nf);

end
% function [Nf,freq,mu]=getFreq(x)
% a=0.01;
% mu=a*max(x);
% N = length(x);
% Nf = 2^nextpow2(N)*16;
% freq = (0:Nf-1)/Nf;
% freq=freq(:);
% end
  
[S_MUSIC1,freqh1,ARspectrum1,mu1]=highResolution(x1);
figure();
plot(freqh1, 10*log10(ARspectrum1));
title('AR Spectral Estimate of signal x1');
xlabel('Normalized Frequency');
ylabel('Power (dB)');
figure();
plot(freqh1, 10*log10(S_MUSIC1));
title('Music Spectral Estimate of signal x1');
xlabel('Normalized Frequency');
ylabel('Power (dB)');
[S_MUSIC2,freqh2,ARspectrum2]=highResolution(x2);
figure();
plot(freqh2, 10*log10(ARspectrum2));
title('AR Spectral Estimate of signal x2');
xlabel('Normalized Frequency');
ylabel('Power (dB)');
figure();
plot(freqh2, 10*log10(S_MUSIC2));
title('Music Spectral Estimate of signal x2');
xlabel('Normalized Frequency');
ylabel('Power (dB)');
[S_MUSIC3,freqh3,ARspectrum3]=highResolution(x3);
figure();
plot(freqh3, 10*log10(ARspectrum3));
title('AR Spectral Estimate of signal x3');
xlabel('Normalized Frequency');
ylabel('Power (dB)');
figure();
plot(freqh3, 10*log10(S_MUSIC3));
title('Music Spectral Estimate of signal x3');
xlabel('Normalized Frequency');
ylabel('Power (dB)');
%%

% [freqh1,mu1]=getFreq(x1);
xS=[x1,x2,x3,x4];
%for x1 signal


Tlab=["x1","x2", "x3", "x4"];

for i=1:length(Tlab)
    xi = xS(:, i);
   [S_MUSIC,freq1,ARspectrum,mu1]=highResolution(xi);
    [u1,axe_freq] = minl1_Fourier(xi,freqh1,mu1);
    figure();
    plot(axe_freq,abs(u1))
    title(['minL1-Fourier Spectrum for ', Tlab(i)]);
    xlabel('Frequency Index');
    ylabel('|u1|');
    grid on;
end    
