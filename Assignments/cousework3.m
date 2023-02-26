N = 8;
k = 1;
fs = 8000;
ts = 1/fs;

p = zeros(1, N-1);

for n = 0:N-1
    t = n *ts;
    x(n+1) = sin(2*pi * 1000 * t) + 0.5*sin(2*pi * 2000 * t+((3*pi)/4));
    %p(n+1) = cos((2 * pi * n * k)/N) - 1i* sin((2 * pi * n * k)/N);
    p(n+1) = exp(-i * pi * 2 * n * k / N);
end

C1 = x * p';
% this matches as per the book (page 66) except the phase is not correct, 90
% vs -90 in book; this is related to the matmul vs element wise and sum


for n = 0:N-1
    t = n *ts;
    x(n+1) = sin(2*pi * 1000 * t) + 0.5*sin(2*pi * 2000 * t+((3*pi)/4));
    for k = 0:N-1
        %p(n+1) = cos((2 * pi * n * k)/N) - 1i* sin((2 * pi * n * k)/N);
        p(k+1, n+1) = exp(-i * pi * 2 * n * k / N);
    end
end

C1 = x * p';
Xk = x * p;


fs = 8000;
ts = 1/fs;
N = 8;
p = zeros(1,N);
m = 1;
for n = 0:(N - 1)
    t = n * ts;
    x = sin(2*pi * 1000 * t) + 0.5*sin(2*pi * 2000 * t+((3*pi)/4));
    p(n+1) = x;
end

for n = 0:N-1
    xn = p(n+1);
    Xm = xn * cos((2 * pi * n * m)/N) - (1i * sin((2 * pi * n * m)/N));
    %Xm = xn * exp(-1i * 2 * pi * n * m /N);
    q(n+1) = Xm;
end





for n = 0:(N - 1)
    t = n * ts;
    x = sin(2*pi * 1000 * t) + 0.5*sin(2*pi * 2000 * t+((3*pi)/4));
    p(n+1) = x;
end

for n = 0:N-1
    for m = 0:N-1
        xn = p(n+1);
        Xm = xn * cos((2 * pi * n * m)/N) - (1i * sin((2 * pi * n * m)/N));
        %Xm = xn * exp(-1i * 2 * pi * n * m /N);
        q(m+1, n+1) = Xm;
    end
end





%{

***********************************************************************
Question 1
***********************************************************************
%}
clear all;

fs = 500;   % sampling rate
ts = 1/fs;  % sample time period
N = 100;
p = zeros(1,N);
t = linspace(0, N*ts, N);
x = sin(2 * pi * 80 * t) + sin(2 * pi * 92 * t);
f_a = zeros(1,N);   % array for the analysis freq.

for n = 0:N-1
    %t = n *ts;
    %x(n+1) = sin(2 * pi * 80 * t) + sin(2 * pi * 92 * t);
    f_a(n+1) = (n * fs)/N;       % analysis freq.
    for k = 0:N-1
        p(k+1, n+1) = exp(-i * pi * 2 * n * k / N);
    end
end

Xk = x * p;


mag = abs(Xk);
phase = angle(Xk);

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
(N/2+1), this is because the second half samples will have a corresponding
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
***********************************************************************
%}
clear all;
close all;

fs = 4;
ts = 1/fs;
N = 0.75*fs;
t = linspace(0,N/fs, N);
x = 5 + 2*cos((2*pi*t) - pi/2) + 3*cos(4*pi*t);

f_a = zeros(1,N);   % array for the analysis freq.

for n = 0:N-1
    %t = n *ts;
    %x(n+1) = sin(2 * pi * 80 * t) + sin(2 * pi * 92 * t);
    f_a(n+1) = (n * fs)/N;       % analysis freq.
    for k = 0:N-1
        p(k+1, n+1) = exp(-i * pi * 2 * n * k / N);
    end
end

Xk = x * p;


mag = abs(Xk);
phase = angle(Xk);

subplot(2,1,1)
stem(f_a, mag)
title('Signal Magnitude Spectra')
grid on
subplot(2,1,2)
stem(f_a,phase)
xlabel('Frequency (Hz)')
ylabel('Phase (Radians)')
title('Signal Phase Spectra')


grid on




