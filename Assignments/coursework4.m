clear all;
close all;

%{
*************************************************************************
Qu1
*************************************************************************

Since this lowpass signal is bandlimited to 40Hz, is satisfies the Nyquist
frequency (fs/2 > B) where fs is the samle rate and B is the highest
frequency band
therefore:
    fs/2 > B    =>     100/2 > 40  =>    50 > 40    =>  True

since the Nyquist equation is satisfied, there will be no alising. However
since this is a band limited signal

Using the equation for the analysis frequency:
   fa = k * fs / N

we get X(50) = 50 * 100 / 400 = 12.5

%}

f0 =40;
fs = 100;
N = 200;
ts = 1/fs;
t = linspace(0, N/ts , N );
x = sin(2*pi*f0*t);

figure(1)
plot(t,x);
ylim([-10,10])

ft = (1/N)* fft(x);

figure(2)
stem(abs(ft))

% add zero padding
x_zeroed = [x, zeros(size(x))];
ft_zeroed = 1/N * fft(x_zeroed);
figure(3)
stem(abs(ft_zeroed))

N_padded = length(x_zeroed);
f_a = zeros(1,N_padded);
for n = 0:N-1
    f_a(n+1) = (n * fs)/N_padded;       % analysis freq.
end

% find analysis frequency for X(50)
X_50 = f_a(51);     % 12.5Hz



%{
*************************************************************************
Qu2

All signals have peaks around 852Hz and 1209Hz, as expected
However the Blackman has performed worst since it has mutiple data points
around the peak, showing signs of leakage. Blackman also has a shallower
trough between the two peaks compared to both rectangular/Hamming - i.e.
hte difference between the peaks and the local minima between the peaks is
less with the Blackman

Rectangular and Hamming though both have more distinct peaks close to 852Hz
and 1209Hz. Either of these would be suitable. Blackman in this case would
not be.

It should be noted that the maximum db magnitudes for all functions are not
exactly 852Hz and 1209Hz, e.g. the value is 844Hz on the first peak. This
could be due to some leakage given that freqencies are close and therefore
hard to distingish.
*************************************************************************
%}

clear all;
close all

% set up the signal and sample time arrays
Fs = 8000 ; 
t = (0:7999)/Fs; 
sig1 = sin(2*pi*852*t);
sig2 = sin(2*pi*1209*t);
x = sig1 + sig2;    % add the signals together

% create the 3 windowing arrays
% Calculate the Fourier transform using rectangular window
X_rect = fft(x, 256);
X_rect_db = db(abs(X_rect));

% fft with Hamming window
w_hamming = hamming(length(x))';
X_hamming = fft(x.*w_hamming, 256);
X_hamming_db = db(abs(X_hamming));

% fft with Blackman window
w_blackman = blackman(length(x))';
X_blackman = fft(x.*w_blackman, 256);
X_blackman_db = db(abs(X_blackman));

% Plot the results
f = Fs*(0:127)/256;
plot(f, X_rect_db(1:128), 'r', f, X_hamming_db(1:128), 'g', f, X_blackman_db(1:128), 'b');
title('FFT: Windowing Comparison');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Rectangular', 'Hamming', 'Blackman');


