close all;
clear all;
clc;
x=load('signals.mat');
xS=[x.x1,x.x2,x.x3,x.x4];
%for x1 signal

Tlab = ["x1", "x2", "x3", "x4"];

for i = 1:length(Tlab)
    [N1, t1, Nf1, freq1, Fx1, Fxw1, xh1, pxx1, f1] = computeNandFft(xS(:,i));

    fig = figure(i);
    
    % --- Time Domain Plot ---
    subplot(3,1,1);
    plot(t1, xS(:,i), "-r", t1, xh1, "b--", 'LineWidth', 2);
    title(['Time domain signal ' Tlab(i)]);
    xlabel('Time');
    ylabel('Amplitude');
    legend({Tlab(i), 'Hamming Windowed'}, 'Location', 'southwest');

    subplot(3,1,2);
    plot(freq1, abs(Fx1), "-r", freq1, abs(Fxw1), "b-", f1, abs(pxx1), "g--", 'LineWidth', 1);
    title(['Periodogram of ' Tlab(i)]);
    xlabel('Frequency');
    ylabel('Power');
    legend({Tlab(i), 'Hamming Windowed', 'Welch'}, 'Location', 'northeast');

    subplot(3,1,3);
    plot(freq1, 10*log10(abs(Fx1) + eps), "-r", ...
         freq1, 10*log10(abs(Fxw1) + eps), "b--", ...
         f1, 10*log10(abs(pxx1) + eps), "g--", 'LineWidth', 1);
    title(['Log Periodogram of ' Tlab(i)]);
    xlabel('Frequency');
    ylabel('Power(dB)');

    filename = ['signal_plot_' char(Tlab(i)) '.pdf'];
    exportgraphics(fig, filename, 'Resolution', 800);
end
    % figure();
    % subplot(3,1,1)
    % plot(t1, xh1);
    % title('Time domain signal '+ Tlab(i));
    % xlabel('Time');
    % ylabel('Amplitude');
    % 
    % subplot(3,1,2);
    % plot(freq1, abs(Fxw1).^2);
    % title('Hamming Window Periodogram of '+ Tlab(i));
    % xlabel('Frequency');
    % ylabel('Power');
    % 
    % subplot(3,1,3);
    % plot(freq1, log(abs(Fxw1).^2));
    % title('log of Hamming Window Periodogram of '+ Tlab(i));
    % xlabel('Frequency');
    % ylabel('Power');
    % 
    

% figure(2);
% subplot(3,1,1)
% plot(t1, xh1);
% title('Time domain signal x1');
% xlabel('Time');
% ylabel('Amplitude');
% 
% subplot(3,1,2);
% plot(freq1, abs(Fxw1).^2);
% title('Hamming Window Periodogram of x1');
% xlabel('Frequency');
% ylabel('Power');
% 
% subplot(3,1,3);
% plot(freq1, log(abs(Fxw1).^2));
% title('log of Hamming Window Periodogram of x1');
% xlabel('Frequency');
% ylabel('Power');
% 
% figure(3);
% subplot(2,1,1)
% plot(f1, pxx1)
% title('Welch periodogram x1');
% xlabel('frequency');
% ylabel('power/frequency');
% 
% subplot(2,1,2);
% plot(f1,10*log10(pxx1));
% title('log of Welch periodogram of x1');
% xlabel('Frequency');
% ylabel('Power/frequency');
%[fig1,fig2]=plotting(x1);
%%
% [N2,t2,Nf2,freq2,Fx2,Fxw2,xh2,pxx2,f2]=computeNandFft(x2);
% figure(4);
% subplot(3,1,1);
% plot(t2, x2);
% title('Time domain signal x2');
% xlabel('Time');
% ylabel('Amplitude');
% 
% subplot(3,1,2);
% plot(freq2, abs(Fx2).^2);
% title('Spactrum of x2');
% xlabel('Frequency');
% ylabel('Magnitude');
% 
% subplot(3,1,3);
% plot(freq2, log(abs(Fx2).^2));
% title('log Spactrum of x2');
% xlabel('Frequency');
% ylabel('Magnitude');
% 
% figure(5);
% subplot(3,1,1);
% plot(t2, xh2);
% title('Time domain signal x2');
% xlabel('Time');
% ylabel('Amplitude');
% subplot(3,1,2);
% plot(freq2, abs(Fxw2).^2);
% title('Hamming Window Periodogram of x2');
% xlabel('Frequency');
% ylabel('Power');
% 
% subplot(3,1,3);
% plot(freq2, log(abs(Fxw2).^2));
% title('log of Hamming Window Periodogram of x2');
% xlabel('Frequency');
% ylabel('Power');
% 
% figure(6);
% subplot(2,1,1)
% plot(f2, pxx2)
% title('Welch periodogram x2');
% xlabel('frequency');
% ylabel('power/frequency');
% 
% subplot(2,1,2);
% plot(f2,10*log10(pxx2));
% title('log of Welch periodogram of x2');
% xlabel('Frequency');
% ylabel('Power/frequency');
% %%
% [N3,t3,Nf3,freq3,Fx3,Fxw3,xh3,pxx3,f3]=computeNandFft(x3);
% figure(7);
% subplot(3,1,1);
% plot(t3, x3);
% title('Time domain signal x3');
% xlabel('Time');
% ylabel('Amplitude');
% 
% subplot(3,1,2);
% plot(freq3, abs(Fx3).^2);
% title('Spactrum of x3');
% xlabel('Frequency');
% ylabel('Magnitude');
% 
% subplot(3,1,3);
% plot(freq3, log(abs(Fx3).^2));
% title('log Spactrum of x3');
% xlabel('Frequency');
% ylabel('Magnitude');
% 
% figure(8);
% subplot(3,1,1);
% plot(t3, xh3);
% title('Time domain signal x3');
% xlabel('Time');
% ylabel('Amplitude');
% 
% subplot(3,1,2);
% plot(freq3, abs(Fxw3).^2);
% title('Hamming Window Periodogram of x3');
% xlabel('Frequency');
% ylabel('Power');
% 
% subplot(3,1,3);
% plot(freq3, log(abs(Fxw3).^2));
% title('log of Hamming Window Periodogram of x3');
% xlabel('Frequency');
% ylabel('Power');
% 
% figure(9);
% subplot(2,1,1)
% plot(f3, pxx3)
% title('Welch periodogram x3');
% xlabel('frequency');
% ylabel('power/frequency');
% 
% subplot(2,1,2);
% plot(f3,10*log10(pxx3));
% title('log of Welch periodogram of x3');
% xlabel('Frequency');
% ylabel('Power/frequency');
% %%
% [N4,t4,Nf4,freq4,Fx4,Fxw4,xh4,pxx4,f4]=computeNandFft(x4);
% figure(10);
% subplot(3,1,1);
% plot(t4, x4);
% title('Time domain signal x4');
% xlabel('Time');
% ylabel('Amplitude');
% 
% subplot(3,1,2);
% plot(freq4, abs(Fx4).^2);
% title('Spactrum of x4');
% xlabel('Frequency');
% ylabel('Magnitude');
% 
% subplot(3,1,3);
% plot(freq4, log(abs(Fx4).^2));
% title('log Spactrum of x4');
% xlabel('Frequency');
% ylabel('Magnitude');
% 
% figure(11);
% subplot(3,1,1);
% plot(t4, xh4);
% title('Time domain signal x4');
% xlabel('Time');
% ylabel('Amplitude');
% 
% subplot(3,1,2);
% plot(freq4, abs(Fxw4).^2);
% title('Hamming Window Periodogram of x4');
% xlabel('Frequency');
% ylabel('Power');
% 
% subplot(3,1,3);
% plot(freq4, log(abs(Fxw4).^2));
% title('log of Hamming Window Periodogram of x4');
% xlabel('Frequency');
% ylabel('Power');
% 
% figure(12);
% subplot(2,1,1)
% plot(f4, pxx4)
% title('Welch periodogram x2');
% xlabel('frequency');
% ylabel('power/frequency');
% 
% subplot(2,1,2);
% plot(f4,10*log10(pxx4));
% title('log of Welch periodogram of x2');
% xlabel('Frequency');
% ylabel('Power/frequency');
% %%
% %figure(1);
% % subplot(3,1,1);
% % plot(t1, x1);
% % title('Time domain signal x1');
% % xlabel('Time');
% % ylabel('Amplitude');
% % 
% % subplot(3,1,2);
% % plot(freq1, abs(Fx1).^2);
% % title('Spactrum of x1');
% % xlabel('Frequency');
% % ylabel('Magnitude');
% % 
% % subplot(3,1,3);
% % plot(freq1, log(abs(Fx1)).^2);
% % title('Spactrum of x1');
% % xlabel('Frequency');
% % ylabel('Magnitude');
% % 
% % figure(2);
% % subplot(3,1,1)
% % plot(t1, xh1);
% title('Time domain signal x1');
% xlabel('Time');
% ylabel('Amplitude');
% subplot(3,1,2);
% plot(freq1, abs(Fxw1).^2);
% title('Hamming Window Periodogram of x1');
% xlabel('Frequency');
% ylabel('Power');
% 
% subplot(3,1,3);
% plot(freq1, log(abs(Fxw1)).^2);
% title('Hamming Window Periodogram of x1');
% xlabel('Frequency');
% ylabel('Power');

% function [N,t,Nf,freq,Fx,Fxw,xh,pxx,f] = computeNandFft(x)
% Fs = 1;
% N = length(x) ;
% t = (0:N-1)/Fs;
% Nf = 2^nextpow2(N)*16;
% freq = 0:Fs/Nf:(Nf-1)/Nf*Fs;
% %fft of signal
% FxP = fft(x,Nf);
% Fx= (1/N*Fs)*abs(FxP).^2;
% %fft of windoed signal
% h = hamming(N);
% xh = x .* h;
% FxwP= fft(xh, Nf);
% Fxw= (1/(N*Fs))*abs(FxwP).^2;
% 
% %welch windowing
% % periodogram_welch = pwelch(x,window,noverlap,Nf,Fs,'twosided');
% % window=hamming(N/4);
% % noverlap=window/4;
% 
% segLen = 2^nextpow2(N / 4); % L
% %segmensize is L if L is small sidelobes will be wider
% %noverlap should be 50% of window
% %nf should be 8*with window size
% window = hamming(segLen);
% % welchWin = 1 - (((0:segLen-1) - (segLen-1)/2) / ((segLen-1)/2)).^2
% noverlap = floor(segLen / 4);
% nfft = 2^nextpow2(segLen)*16;   % bigger for better freq resolution
% 
% [pxx, f] = pwelch(x, window,noverlap, nfft, Fs,'twosided');
% %Fx_welch = fft(pxx, Nf);