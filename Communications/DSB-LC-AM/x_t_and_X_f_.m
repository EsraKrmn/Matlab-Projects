f = 4000;             % defining frequency
t = 0:(1/f):1-(1/f);  % defining time scale

% message signal m(t): 
mt = (5*cos(2*pi*100*t)) + (10*cos(2*pi*200*t));
xt = mt + 2*cos(2*pi*1000*t);  % x(t) function

% Plot x(t)
subplot(2, 1, 1); % Create axes for graphs and put xt to 1th
plot(t,xt);
title('x(t) signal');
xlabel('t(s)');
ylabel('x(t)');

% Find Fast fourier transform of x(t)
fourier_xt = fft(xt)/length(xt);  
freq = -f/2:(f/length(xt)):f/2-(f/length(xt));
Xf = abs(fftshift(fourier_xt));

% Plot FFT of x(t)
subplot(2, 1, 2);     % Create axes for graphs and put Xf to 2nd
stem(freq, abs(Xf));  % stem for plotting discrete sequence data 
title('FFT of x(t) : X(f)');
xlabel('f(Hz)');
ylabel('|X(f)|');