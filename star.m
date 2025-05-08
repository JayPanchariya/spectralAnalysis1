close all;
clear all;
clc;
load ('stardata_2022.mat');

N=length(x);
Nf=4096;
% freq = 0:Fs/Nf:(Nf-1)/Nf*Fs;
xw=x.*win;

[freq,yfft] = compute_fft_and_shift(xw,Nf,Fs);
xf = (1 / (Fs * N)) * abs(yfft).^2;

% single-sided and adjust power
yfft = yfft(Nf/2+1 : Nf);
xf1  = xf(Nf/2+1 : Nf);
xf1(2:end-1) = 2 * xf1(2:end-1);  
fq   = freq(Nf/2+1 : Nf);




fig1=figure;
subplot(3,1,1)
plot(t,xw)
title('Star Data');
xlabel('Time');
ylabel('Signal')

subplot(3,1,2);
plot(fq,abs(yfft));
title('Spactrum of Star Data');
xlabel('Frequency');
ylabel('Magnitude');
xlim([5 25])

subplot(3,1,3);
plot(fq,10*log(abs(yfft)));
title('Log Spactrum of Star Data');
xlabel('Frequency');
ylabel('Magnitude');
xlim([5 25])
ylim([-50 50])
exportgraphics(fig1, 'stardataspectrum1.pdf', 'Resolution', 800);

fig2=figure;
subplot(3,1,1)
plot(t,x)
xlabel('Time');
ylabel('Signal')
title('Star Data');

subplot(3,1,2);
plot(fq,(xf1));
title('Periodogram of  stardata');
xlabel('Frequency');
ylabel('Power');
xlim([5 25])

subplot(3,1,3);
plot(fq,10*log(abs(xf1)));
title('Log Periodogram of  stardata');
xlabel('Frequency');
ylabel('Power');
xlim([5 25])
ylim([-200 10])

exportgraphics (fig2 ,'stardataPeriodogram1.pdf' , 'Resolution', 800);


yP = yfft(length(yfft)/2: length(yfft));
fp = freq(length(yfft)/2: length(yfft));

ypo=yP;

fpo=fp;
%%

residual = xw;
threshold = 0.0001;
amp = Inf;  
iter=0;

while amp > threshold
    [freq, yfft] = compute_fft_and_shift(residual, Nf, Fs);
    xf = (1 / (Fs * Nf)) * abs(yfft).^2; 

 
    center = floor(length(yfft)/2);
    yP = xf(center+1:end);   % Positive frequencies
    fp = freq(center+1:end); % Corresponding frequencies
    length(residual)
    [~, idx] = max(yP);
   
    l=length(yfft);
    amp = abs(yP(idx));
    fi = fp(idx);  


    [alpha, beta] = estim_cos_sin(residual, t, fi);
    current = alpha * cos(2 * pi * fi * t) + beta * sin(2 * pi * fi * t);

    residual = (residual-current').*win;
    % Plot spectrum
    fig3 = figure;
    plot(fp, yP, '-', 'Marker', 'o', 'MarkerIndices', idx,'MarkerEdgeColor', 'r', 'LineWidth', 1, 'Color', 'b');
    title(['Signal Spectrum at iteration ' num2str(iter)]);
    xlabel('Frequency');
    ylabel('Magnitude');
    spectrum_filename = sprintf('iterSpectrum_%d.pdf', iter);
    exportgraphics(fig3, spectrum_filename, 'Resolution', 800);
    close(fig3);  

    % Plot residual
    fig4 = figure;
    plot(t, xw, '-r', t, residual, '--b');
    title(['Residual at iteration ' num2str(iter)]);
    xlabel('Time');
    ylabel('Signal');
    legend({'Signal', 'Residual Signal'}, 'Location', 'northeast');
    residual_filename = sprintf('iterResidual_%d.pdf', iter);
    exportgraphics(fig4, residual_filename, 'Resolution', 800);
    close(fig4);
    
    iter=iter+1;
end

sigN=xw-(residual);
[freq1,cleanSig] = compute_fft_and_shift(sigN,Nf,Fs);
yP = cleanSig(length(cleanSig)/2: length(yfft));
resP = (1 / (Fs * N)) * abs(cleanSig).^2;

% Keep single-sided and adjust power
resP= resP(Nf/2+1 : Nf);
resP(2:end-1) = 2 * resP(2:end-1); 
freq1 = freq1(Nf/2+1 : Nf);
fig5=figure;


subplot(3,1,1)
plot(t,xw,'-r',t,sigN,'--b',LineWidth=1)
xlabel('Time');
ylabel('Signal')
xlim([0 5])
legend({'Signal', 'Noiseless Signal'},'Location','northeast')

subplot(3,1,2);
plot(fq,abs(xf1),'-r',freq1, abs(resP),'--b', LineWidth=1);
title('Periodogram of stardata');
xlabel('Frequency');
ylabel('Power');
xlim([5 25])


subplot(3,1,3);
plot(fq,10*log(abs(xf1)),'-r',freq1,10*log(abs(resP)),'--b',LineWidth=1);
title('Log Periodogram of stardata');
xlabel('Frequency');
ylabel('Power');
xlim([5 25])
ylim([-200 10])

exportgraphics (fig5 ,'NoiselessSig.pdf' , 'Resolution', 800);

%%
% figure();
% plot(freq1,10*log(abs(yfft)),'-r',freq1,10*log(abs(cleanSig)),'--b',LineWidth=2)
% legend({'Signal Spactrum', 'Noiseless Signal spactrum'},'Location','northeast')

% subplot(3,2,2)
% plot(t,sigN)
% xlabel('Time');
% ylabel('Noiseless Signal')
% 
% subplot(3,2,4);
% plot(freq1,abs(cleanSig));
% title('Spactrum of Noiseless Signal');
% xlabel('Frequency');
% ylabel('Magnitude');
% 
% subplot(3,2,6);
% plot(freq1,10*log(abs(cleanSig)));
% title('Log Spactrum of Noiseless Signal');
% xlabel('Frequency');
% ylabel('Magnitude');
% fig3=figure;
    % plot(fp, yP, '-', 'Marker', 'o', 'MarkerIndices', idx, 'MarkerEdgeColor', ...
    % 'r','LineWidth', 1, 'Color', 'b');
    % title(['Signal Spactrum at iteration ' num2str(iter)]);
    % xlabel('Frequency');
    % ylabel('Power')
    % exportgraphics (fig3 , 'iterSpactrum.pdf' , 'Resolution', 800);
    % 
    % 
    % drawnow;
    % 
    % fig4=figure;
    % plot(t,xw,'-r',t,residual,'--b' )
    % title(['Residual at ' num2str(iter)] )
    % xlabel('Time');
    % ylabel('Signal')
    % legend({'Signal', 'Residual Signal'},'Location','northeast')
    % exportgraphics (fig4 ,'iterResSig.pdf' , 'Resolution', 800);