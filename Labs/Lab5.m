% Windowing
%Week 5 Lab work

%% Ex1
Fs = 512 ; 
t = (0:511)/Fs; 
f = (0:2048)/2048*(Fs/2);
recwin=rectwin(512);
fftrectwin =fft(recwin, 4096) ;
plot(f, db(abs(fftrectwin(1:2049)))); % db converts to decibles
xlabel ('Frequency(Hz)'); % labels x-axis of spectral plot
ylabel ('Magnitude Spectrum (dB)'); % labels y-axis of spectral plot
% i) ~13
% ii) ~ 2

%% Ex2
hamwin=hamming(512); % Generates a Hamming window with N= 512 samples
ffthamwin = fft(hamwin, 4096) ;
figure; plot(t, hamwin)
figure ; plot(f, db(abs(ffthamwin(1:2049))));

blackwin=blackman(512); % Generates a Blackman window with N= 512 samples
fftblackwin = fft (blackwin, 4096);
figure; plot(t, blackwin)
figure ; plot(f, db(abs(fftblackwin(1:2049))));

%% Ex3
close all
S1=sin(2*pi*20*t);
sig=S1;
sig_rectwin=recwin'.*sig;
fftrectsig=fft(sig_rectwin, 512);
f = (0:256)/256*(Fs/2);
plot(f, db(abs(fftrectsig(1:257))))
xlabel ('Frequency(Hz)');
ylabel ('Magnitude Spectrum (dB)');
%  48.1648
S2=sin(2*pi*100.5*t);
sig=S1+S2;
sig_rectwin=recwin'.*sig;
fftrectsig=fft(sig_rectwin, 512);
f = (0:256)/256*(Fs/2);
figure; plot(f, db(abs(fftrectsig(1:257))))
xlabel ('Frequency(Hz)');
ylabel ('Magnitude Spectrum (dB) - rectangular');
% 48.1649

%%
sig_hamwin = hamwin'.*sig;
ffthamsig = fft (sig_hamwin , 512);
figure; plot(f, db(abs(ffthamsig(1:257))))
xlabel ('Frequency(Hz)');
ylabel ('Magnitude Spectrum (dB) - hamming');


%%
sig_black = blackwin'.*sig;
fftblacksig = fft (sig_black , 512);
figure; plot(f, db(abs(fftblacksig(1:257))))
xlabel ('Frequency(Hz)');
ylabel ('Magnitude Spectrum (dB) - blackman');


%%Q4
S3=0.05*sin(2*pi*17*t);
sig=S1+S2+S3;
%Apply rectangular window

sig_rectwin=recwin'.*sig;
fftrectwin =fft(sig_rectwin, 512) ;
figure;
plot(f, db(abs(fftrectwin(1:257)))); % db converts to decibles
xlabel ('Frequency(Hz)');
ylabel ('Magnitude Spectrum (dB) - rectangular');
title("rectangular window")

%Blackman window
blackwin=blackman(512); % Generates a Blackman window with N= 512 samples
sig_black = blackwin.*sig';
fftblackwin = fft (sig_black, 512);
figure; plot(t, blackwin)
title("Blackwin")
figure ; plot(f, db(abs(fftblackwin(1:257))));
xlabel ('Frequency(Hz)');
ylabel ('Magnitude Spectrum (dB) - blackman');

%Hamming window
hamwin=hamming(512); % Generates a Hamming window with N= 512 samples
sig_hamwin = hamwin.*sig';
ffthamwin = fft(sig_hamwin, 512) ;
figure; plot(t, hamwin)
title("rectangular window")
figure ; plot(f, db(abs(ffthamwin(1:257))));
xlabel ('Frequency(Hz)');
ylabel ('Magnitude Spectrum (dB) - hamming');