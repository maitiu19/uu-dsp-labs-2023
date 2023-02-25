N = 8;
k = 1;
fs = 8000;
ts = 1/fs;

p = zeros(1, N-1);

for n = 0:N-1
    t = n *ts;
    x(n+1) = sin(2*pi * 1000 * t) + 0.5*sin(2*pi * 2000 * t+((3*pi)/4));
    p(n+1) = cos((2 * pi * n * k)/N) - 1i* sin((2 * pi * n * k)/N);
    %p(n+1) = exp(-i * pi * 2 * n * k / N);
end

C1 = x * p';



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






clear all;

fs = 500;
ts = 1/500;
N = 100;
p = zeros(1,N);
for n = 0:(N - 1)
    t = n * ts;
    x = sin(2 * pi * 80 * t) + sin(2 * pi * 92 * t);
    p(n+1) = x;
end

for n = 0:(N-1)
    m = n;
    xn = p(n+1);
    Xm = xn * cos((2 * pi * n * m)/N) - (1i * sin((2 * pi * n * m)/N));
    q(n+1) = Xm;
end
