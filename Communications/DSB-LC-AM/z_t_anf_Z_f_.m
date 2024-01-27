f = 4000;             % defining frequency
t = 0:(1/f):1-(1/f);  % defining time scale

% message signal m(t): 
mt = (5*cos(2*pi*100*t)) + (10*cos(2*pi*200*t));
xt = mt + 2*cos(2*pi*1000*t);  % x(t) function
zt = (60*xt) + (xt.^2);        % z(t) function

% Plot z(t)
subplot(2, 1, 1); % Create axes for graphs and put zt to 1th
plot(t,zt);
title('z(t) signal');
xlabel('t(s)');
ylabel('z(t)');

% Find Fast fourier transform of z(t)
fourier_zt = fft(zt)/length(zt);  
freq = -f/2:(f/length(zt)):f/2-(f/length(zt));
Zf = abs(fftshift(fourier_zt));

% Plot FFT of z(t)
subplot(2, 1, 2);     % Create axes for graphs and put Zf to 2nd
stem(freq, abs(Zf));  % stem for plotting discrete sequence data 
title('FFT of z(t) : Z(f)');
xlabel('f(Hz)');
ylabel('|Z(f)|');