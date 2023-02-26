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

Due to DFT symatry
***********************************************************************
%}


