f = 600;              % defining frequency
t = 0:(1/f):1-(1/f);  % defining time scale
%t = 0:0.0001:0.05;   %It is for seeing mt clearly

% message signal m(t): 
mt = (5*cos(2*pi*100*t)) + (10*cos(2*pi*200*t));

% Plot message signal
subplot(2, 1, 1); % Create axes for graphs and put mt to 1th
plot(t,mt);  
title('Message signal');
xlabel('t(s)');
ylabel('m(t)');

% Find Fast fourier transform of m(t)
Mf = fft(mt)/length(mt);
freq = linspace(-f/2, f/2, length(Mf));

% Plot FFT of Message signal
subplot(2, 1, 2);  % Create axes for graphs and put Mf to 2nd
stem(freq, abs(Mf));  % stem for plotting discrete sequence data 
title('FFT of m(t) : M(f)');
xlabel('f(Hz)');
ylabel('|M(f)|');
