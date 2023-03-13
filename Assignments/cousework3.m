%{
***********************************************************************
Question 1

Comments on output:
Shape of the Amplitude sketch is in line with the theory, of the 100 points
sampled, the plot is symetrical at the halfway point - this is because for
any X(k) where 0 <= k <= (N/2 -1), X(k) is the complex conjugate of X(N-k),
so the magnitudes will be equal, hence the symmetry

Regarding the sample peaks - the maximum |X(k)| is close to 50 - the reason
for this is that when a signal has a frequency of f<2*Fs (which is the
case here in order to satisty the Nyquist freq.), the max. magnitude of the
DFT component will be M = A * N/2 where A is the signal amplitude.

Finally, the frequency bins, the peak magnitude of ~50 at f = 80Hz (see
plot), does not mean that the highest magnitude component in X(k) has a
frequency of 80Hz. Since the x-axis is the analysis frequency (Fa), given by
equation:
Fa(k) = k*fs / N => Fa(18) = 18 * 500 / 100 = 90Hz

Matlab logic for this is below the plot, 
I'm not 100% certain it  is correct, there might be a mistake with the
matlab index conversion to k, but the point is that the frequency of the
largest magnitude component will depend on the sample frequency, Fs

Above theory from:
R. G. Lyons, “3.5 DFT Frequency Axis,” in Understanding Digital Signal
 Processing, Upper Saddle River etc.: Pearson Education International,
 2013, pp. 74–77. 
***********************************************************************
%}

clear all;
close all;

fs = 500;   % sampling rate
ts = 1/fs;  % sample time period
N = 100;    % number of samples
p = zeros(1,N); % define array for the exp. component of Xk(n)
t = linspace(0, N*ts, N);   % time array of N samples
x = sin(2 * pi * 80 * t) + sin(2 * pi * 92 * t); %  x(n) component of Xk(n)
f_a = zeros(1,N);   % array for the analysis freq.

for n = 0:N-1
    f_a(n+1) = (n * fs)/N;       % analysis freq.
    for k = 0:N-1
        p(k+1, n+1) = exp(-i * pi * 2 * n * k / N);
    end
end

Xk = x * p;

mag = abs(Xk);
phase = angle(Xk);

%{
% extra code to confirm frequency bin theory:
    [xk_max, xk_max_index] = max(mag);

result is xk_max_index = 85, but k = 0 at mag(1), so
since X(k) = |X(N - k)| the max will also be found at 100 - 85 = 15
this is 85 in matlab, but k = 84 and evaluates to:
k = |N+1 - k| = mag(17) = 49.6438, which corresponds to k(18)

%}

% plot magnitude spectra
figure(1)
stem(f_a, mag)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Discrete Fourier Transform of Signal')
grid on


%{

***********************************************************************
Question 2

Due to DFT symmmetry, we only need the first half of samples from a signal
(N/2 - 1), this is because the second half samples will have a corresponding
conjugate in the first half of samples, i.e.:
X(m) =  X^*(N - m) , where X* is the conjugate, in other words:

x = exp(-j * 2pi) ==> x* = exp(j * 2pi)

therefore the sample:
X(0) =  5
X(1) =  2 − j5
X(2) =  −11.8 + j1.8
X(3) =  12.85 + j1.2
X(4) =  −1 − 3.4j
X(5) =  0.5 − 0.866j
X(6) =  6 − 1.9j
X(7) =  12.8 + 5j

will continue as:
X(8) =  5
X(9) =  2 + j5
X(10) =  −11.8 - j1.8
X(11) =  12.85 - j1.2
X(12) =  −1 + 3.4j
X(13) =  0.5 + 0.866j
X(14) =  6 + 1.9j
X(15) =  12.8 - 5j
***********************************************************************
%}


%{

***********************************************************************
Question 3

Comments on the output:
The sample rate is 4Hz, i.e. 4 samples per second. 
The lowest signal frequency is 1Hz fom the first term:
    cos(2*pi*t - 90deg). 
So we should have N samples given by:
    N = freq_min * fs = 1 * 4 --> 4 samples minimum. 

But we're asked to sample over 3/4s at a sample rate of 4/s, ie. 3 samples
are taken. This isn't enough and results in signal leakage (e.g. at bin 0,
we should see the DC componant only, magnitude 5. Instead it is 4.8,
perhaps due to leakage of the other signals which have lower amplitudes.

We should have at least 4 samples here, so a longer time period would be
needed if using fs = 4Hz. 
Ideally, increasing the sample rate also, since fs = 4Hz just about
satisfies the Nyquist frequency:
    fs => 2 * f_max = 4Hz
(f_max, highest frequency = 2Hz from the second cosine component)
***********************************************************************
%}

clear all;
close all;

fs = 4;
N = 0.75*fs;
t = linspace(0,N/fs, N);
x = 5 + 2*cos((2*pi*t) - pi/2) + 3*cos(4*pi*t);

f_a = zeros(1,N);   % array for the analysis freq.

for n = 0:N-1
    f_a(n+1) = (n * fs)/N;       % analysis freq.
    for k = 0:N-1
        p(k+1, n+1) = exp(-i * pi * 2 * n * k / N);
    end
end

Xk = (1/ N) * x * p;

% magnitude is normalised
mag = abs(Xk);
phase = angle(Xk);

% plot the magnitude and phase spectra
subplot(2,1,1)
stem(f_a, mag)
title('Signal Magnitude Spectra')
ylabel('Amplitude')
grid on
subplot(2,1,2)
stem(f_a,phase)
xlabel('Frequency (Hz)')
ylabel('Phase (Radians)')
title('Signal Phase Spectra')
grid on



