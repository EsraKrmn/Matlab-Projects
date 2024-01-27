f = 4000;              % defining frequency
t = 0:(1/f):1-(1/f);  % defining time scale

% message signal m(t): 
mt = (5*cos(2*pi*100*t)) + (10*cos(2*pi*200*t));
xt = mt + 2*cos(2*pi*1000*t);  % x(t) function
zt = (60*xt) + (xt.^2);        % z(t) function
freq = -f/2:(f/length(zt)):f/2-(f/length(zt));

BW = 400;    % Bandwidth
fc = 1000;   % c enter frequency
Gain = 1;    % Gain
filter_freq = abs(freq);
BGF_down = (filter_freq >= fc-(BW/2)); % Defining lower limit
BGF_up = (filter_freq <= fc+(BW));     % Defining upper limit
BGF = (BGF_down & BGF_up) * 1;   % Combine each of them and multiply Gain=1

% Apply bandpass filter to z(t) to get y(t) and Y(f)
Yf = (fft(zt)/length(zt)) .* fftshift(BGF);
yt = ifft(Yf);

%plot y(t)
subplot(2, 1, 1);
plot(t, yt);
title('Output Signal y(t)');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot FFT of y(t)
subplot(2, 1, 2);
stem(freq, Yf);
title('FFT of y(t) : Y(f)');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');

