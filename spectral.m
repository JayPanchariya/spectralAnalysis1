close all;
clear all;
clc;
x=load('signals.mat');
xS=[x.x1,x.x2,x.x3,x.x4];


function [N, t, Nf, freq, Fx, Fxw, xh, pxx, f] = computeNandFft(x)
    Fs = 1;                       
    N = length(x);              
    t = (0:N-1) / Fs;            

    % Frequency resolution with zero-padding (high-res FFT)
    Nf = 2^nextpow2(N) * 16;
    f_full = (0:Nf-1) * Fs / Nf; % Full FFT frequency vector

    % === Raw FFT Power Spectrum ===
    FxP = fft(x, Nf);
    Fx = (1 / (Fs * N)) * abs(FxP).^2;  % Power spectral density

    % Single-sided spectrum (0 to Nyquist)
    Fx = Fx(1:Nf/2+1);
    Fx(2:end-1) = 2 * Fx(2:end-1);      % Double for power conservation
    freq = f_full(1:Nf/2+1);           % Match single-sided frequencies

    % === Windowed FFT (Hamming) ===
    h = hamming(N);
    xh = x .* h;                        % Windowed signal
    U = sum(h.^2) / N;                 % Normalization factor

    FxwP = fft(xh, Nf);
    Fxw = (1 / (Fs * N * U)) * abs(FxwP).^2;

    Fxw = Fxw(1:Nf/2+1);
    Fxw(2:end-1) = 2 * Fxw(2:end-1);   % Single-sided adjustment

    % === Welch Periodogram ===
    segLen = 2^nextpow2(N / 4);       % Segment length for Welch
    window = hamming(segLen);        
    noverlap = floor(segLen / 4);     
    nfft = 2^nextpow2(segLen) * 16;   % High resolution

    [pxx, f] = pwelch(x, window, noverlap, nfft, Fs, 'psd');

end



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

    % --- Periodogram (Linear) ---
    subplot(3,1,2);
    plot(freq1, abs(Fx1), "-r", freq1, abs(Fxw1), "b-", f1, abs(pxx1), "g--", 'LineWidth', 1);
    title(['Periodogram of ' Tlab(i)]);
    xlabel('Frequency');
    ylabel('Power');
    legend({Tlab(i), 'Hamming Windowed', 'Welch'}, 'Location', 'northeast');

    % --- Log Periodogram (dB) ---
    subplot(3,1,3);
    plot(freq1, 10*log10(abs(Fx1) + eps), "-r", ...
         freq1, 10*log10(abs(Fxw1) + eps), "b--", ...
         f1, 10*log10(abs(pxx1) + eps), "g--", 'LineWidth', 1);
    title(['Log Periodogram of ' Tlab(i)]);
    xlabel('Frequency');
    ylabel('Power(dB)');

    % --- Save with Unique Filename ---
    filename = ['signal_plot_' char(Tlab(i)) '.pdf'];
    exportgraphics(fig, filename, 'Resolution', 800);
end

%%
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
% ARspectrum= ARspectrum(Nf/2+1 : Nf);
% ARspectrum(2:end-1) = 2 * ARspectrum(2:end-1); 
% S_MUSIC  = S_MUSIC(Nf/2+1 : Nf);
% S_MUSIC(2:end-1) = 2 * S_MUSIC(2:end-1);  % Double all but DC and Nyquist
% freq   = freq(Nf/2+1 : Nf);
end

for i = 1:length(Tlab)
    xi = xS(:,i);
    
    % Compute classical spectra
    [N1, t1, Nf1, freq1, Fx1, Fxw1, xh1, pxx1, f1] = computeNandFft(xi);
    
    % Compute AR and MUSIC
    [S_MUSIC, freq_hr, ARspectrum, mu] = highResolution(xi);
    halfIdx = 1:Nf1/2+1;
    ARspectrum = ARspectrum(halfIdx);
    S_MUSIC = S_MUSIC(halfIdx);
    freq_hr = freq_hr(halfIdx);

    % Normalize
    Fx1 = Fx1 / max(Fx1);
    Fxw1 = Fxw1 / max(Fxw1);
    pxx1 = pxx1 / max(pxx1);
    ARspectrum = ARspectrum / max(ARspectrum);
    S_MUSIC = S_MUSIC / max(S_MUSIC);

    % --- High-Resolution Spectrum Plot ---
    fig_hr = figure(100 + i); clf;
    plot(freq_hr, ARspectrum, 'k', 'LineWidth', 1.5); hold on;
    plot(freq_hr, S_MUSIC, 'm--', 'LineWidth', 1.5);
    title(['High Resolution Spectrum for ' Tlab(i)]);
    xlabel('Frequency'); ylabel('Power / Pseudo-Power');
    legend({'AR Spectrum', 'MUSIC Spectrum'}, 'Location', 'northeast');
    grid on; 
    xlim([0 0.5]);

    % Save high-resolution spectrum plot
    filename_hr = ['HR_spectrum_' char(Tlab(i)) '.pdf'];
    exportgraphics(fig_hr, filename_hr, 'Resolution', 800);

    % --- Full Spectrum Comparison Plot ---
    fig = figure(200 + i); clf;
    plot(freq1, Fx1, '-r', freq1, Fxw1, '-b',1, pxx1, 'g--', freq_hr, ARspectrum, 'k:',freq_hr, S_MUSIC, 'm--' ,'LineWidth', 1);
    title(['Normalized Spectrum Comparison - ' Tlab(i)]);
    xlabel('Frequency'); ylabel('Normalized Power');
    legend({'Raw FFT', 'Hamming', 'Welch', 'AR', 'MUSIC'}, 'Location', 'northeast');
    grid on;
    xlim([0 0.5]);

    % Save combined spectrum comparison plot
    filename_combined = ['combined_spectrum_' char(Tlab(i)) '.pdf'];
    exportgraphics(fig, filename_combined, 'Resolution', 800);
end

for i = 1:length(Tlab)
    xi = xS(:, i);

    [S_MUSIC, freq1, ARspectrum, mu1] = highResolution(xi);

    
    [u1, axe_freq] = minl1_Fourier(xi, freq1, mu1);

   
    fig_minl1 = figure(300 + i); clf;
    plot(axe_freq, abs(u1), 'b', 'LineWidth', 1.5);
    title(['minL1-Fourier Spectrum for ', Tlab(i)]);
    xlabel('Frequency Index');
    ylabel('|u1|');
    grid on;

    
    filename_minl1 = ['minL1_spectrum_' char(Tlab(i)) '.pdf'];
    exportgraphics(fig_minl1, filename_minl1, 'Resolution', 800);
end  

% for i = 1:length(Tlab)
%     xi = xS(:,i);
% 
%     % --- FFTs and Welch ---
%     [N1, t1, Nf1, freq1, Fx1, Fxw1, xh1, pxx1, f1] = computeNandFft(xi);
% 
%     % --- AR and MUSIC ---
%     [S_MUSIC, freq_hr, ARspectrum, mu] = highResolution(xi);
% 
%     % Keep single-sided spectra
%     halfIdx = 1:Nf1/2+1;
%     ARspectrum = ARspectrum(halfIdx);
%     S_MUSIC = S_MUSIC(halfIdx);
%     freq_hr = freq_hr(halfIdx);
% 
%     % === Plot all in one figure ===
%     fig = figure(i);
% 
%     % --- Time Domain Plot ---
%     subplot(3,1,1);
%     plot(t1, xi, "-r", t1, xh1, "b--", 'LineWidth', 2);
%     title(['Time Domain Signal: ' Tlab(i)]);
%     xlabel('Time'); ylabel('Amplitude');
%     legend({Tlab(i), 'Hamming Windowed'}, 'Location', 'southwest');
% 
%     % --- Linear Power Spectrum ---
%     subplot(3,1,2);
%     plot(freq1, abs(Fx1), "-r", ...
%          freq1, abs(Fxw1), "b-", ...
%          f1, abs(pxx1), "g--", ...
%          freq_hr, abs(ARspectrum), "k:", ...
%          freq_hr, abs(S_MUSIC), "m--", ...
%          'LineWidth', 1);
%     title(['Power Spectrum Comparison (Linear) - ' Tlab(i)]);
%     xlabel('Frequency'); ylabel('Power');
%     legend({'Raw FFT', 'Hamming', 'Welch', 'AR', 'MUSIC'}, 'Location', 'northeast');
% 
%     % --- Log Power Spectrum ---
%     subplot(3,1,3);
%     plot(freq1, 10*log10(abs(Fx1) + eps), "-r", ...
%          freq1, 10*log10(abs(Fxw1) + eps), "b--", ...
%          f1, 10*log10(abs(pxx1) + eps), "g--", ...
%          freq_hr, 10*log10(abs(ARspectrum) + eps), "k:", ...
%          freq_hr, 10*log10(abs(S_MUSIC) + eps), "m--", ...
%          'LineWidth', 1);
%     title(['Power Spectrum Comparison (dB) - ' Tlab(i)]);
%     xlabel('Frequency'); ylabel('Power (dB)');
%     legend({'Raw FFT', 'Hamming', 'Welch', 'AR', 'MUSIC'}, 'Location', 'northeast');
% 
%     % --- Optional: Save PDF ---
%     filename = ['spectrum_comparison_' char(Tlab(i)) '.pdf'];
%     exportgraphics(fig, filename, 'Resolution', 800);
% end



%for x1 signal
% 
% function [N,t,Nf,freq,Fx,Fxw,xh,pxx,f] = computeNandFft(x)
% Fs = 1;
% N = length(x);
% t = (0:N-1) / Fs;
% 
% %Frequency resolution 
% Nf = 2^nextpow2(N) * 16;
% freq = (0:Nf-1) * Fs / Nf;
% 
% 
% % FFT and single sided PSD
% FxP = fft(x, Nf);
% Fx = (1 / (Fs * N)) * abs(FxP).^2;
% 
% % Keep single-sided and adjust power
% Fx = Fx(Nf/2+1 : Nf);
% Fx(2:end-1) = 2 * Fx(2:end-1);  % Double non-DC/Nyquist
% freq = freq(Nf/2+1 : Nf);
% 
% % Windowed FFT (Hamming)
% h = hamming(N);
% xh = x .* h;
% U = sum(h.^2) / N;  % Normalization factor to compensate window power loss
% 
% FxwP = fft(xh, Nf);
% Fxw = (1 / (Fs * N * U)) * abs(FxwP).^2;  % PSD scaling with window correction
% 
% % Keep single-sided spectrum
% Fxw = Fxw(Nf/2+1 : Nf);
% Fxw(2:end-1) = 2 * Fxw(2:end-1);  
% freq = freq(Nf/2+1 : Nf);  
% 
% % Welch Periodogram 
% segLen = 2^nextpow2(N / 4);      % Segment length L
% window = hamming(segLen);        % Window for Welch
% noverlap = floor(segLen / 4);    % 25% overlap
% nfft = 2^nextpow2(segLen) * 16;  % More zero padding for resolution
% 
% [pxx, f] = pwelch(x, window, noverlap, nfft, Fs, 'psd');  % Proper PSD output
% end
  % 
  % === Zoom-in Inset (subplot 2) ===
  %   axZoom = subplot(2,1,2);
  %   zoom_range = [0.05 0.2];  % adjust based on signals
  %   zoom_idx = find(freq1 >= zoom_range(1) & freq1 <= zoom_range(2));
  %   plot(freq1(zoom_idx), Fx1(zoom_idx), '-r', ...
  %        freq1(zoom_idx), Fxw1(zoom_idx), '-b', ...
  %        f1, pxx1, 'g--', %...
  %        freq_hr, ARspectrum, 'k:', ...
  %        freq_hr, S_MUSIC, 'm--', ...
  %        'LineWidth', 1);
  %   title('Zoom on Main Peak Region');
  %   xlabel('Frequency'); ylabel('Normalized Power');
  %   legend({'Raw FFT', 'Hamming', 'Welch', 'AR', 'MUSIC'}, 'Location', 'northeast');
  %   grid on;
  %   xlim(zoom_range);
  % 
  %   % Save high-quality PDF
  %   filename = ['spectrum_zoom_comparison_' char(Tlab(i)) '.pdf'];
  %   exportgraphics(fig, filename, 'Resolution', 800);