% Qu1

x = [1, 1, 0, 0, 0];
N = 5;
k = 1;

p = zeros(1, N-1);

for n = 0:N-1
    p(n+1) = exp(-i * pi * 2 * n * k / N);
end

C1 = 1/N * x * p';
amp_1 = abs(C1);
phase_1 = angle(C1);

disp(amp_1);
disp(phase_1);

clear all;
close all;

% adapting for all coefficents
x = [1, 1, 0, 0, 0];
N = 5;

Xk_n = zeros(1, N-1);
p = zeros(1, N-1);

for n = 0:N-1
    for k = 0:N-1
        p(k+1, n+1) = exp(-1i * pi * 2 * n * k / N);
    end
end

Xk = 1/N * x * p;
amp_1 = abs(Xk(1));
phase_1 = angle((Xk(1)));
   
disp(amp_1);
disp(phase_1);
   
clear all;

% Qu3
x = [1, 2, 1, 0, -1, 1];
N = length(x);
fs = 1000; % 1kHz sample rate
delta_fs = 1000/N;
k = 2;
f = k * delta_fs;

for n = 0:N-1
    p(n+1) = exp( -1i * 2 * pi * n * k /N);
end

Xk_2 = x * p';

amp_2 = abs(Xk_2);
phase_2 = angle((Xk_2));
   
disp(amp_2);
disp(phase_2);

clear all;
% Qu 4
x = [1,2,1,0,-1,1];
N = length(x);
fs = 1000; % 1kHz sample rate
delta_fs = 1000/N;
for n = 0:N-1
    for k = 0:N-1
        p(k+1, n+1) = exp( -1i * 2 * pi * n * k /N);
        f(k+1) = (k+1) * delta_fs;
    end
end

Xk_n = x * p;

amp = abs(Xk_n);
phase = angle((Xk_n));
   
disp(amp);
disp(phase);

subplot(2,1,1)
stem(f, amp)
subplot(2,1,2)
stem(f,phase)



