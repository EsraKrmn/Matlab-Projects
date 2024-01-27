fs = 100000;   % Sampling frequency
Ts = 1/fs ;    % Sampling period
t = -2:Ts:6 ;  % Time vector

% Define the signal m(t)
m_t = rectpuls((t-1)/2) - rectpuls((t-3)/2);

% Plot m(t) in the time domain
figure;
subplot(2, 1, 1);
plot(t, m_t);
title('Message Signal m(t)');
xlabel('t');
ylabel('m(t)');
grid on;

% Calculate and plot the magnitude spectrum of m(t)
n = length(m_t);               % For normalize 
fshift = (-n/2:n/2-1)*(fs/n);  % Calculate f parameter
M_f = abs(fftshift(fft(m_t)))*Ts;  % Calculate M(f)

subplot(2, 1, 2);
plot(fshift, M_f);
title('Magnitude Spectrum of m(t)');
xlabel('Frequency (Hz)');
ylabel('|M(f)|');
grid on;

% Change the x-axis limits
xlim([-50, 50]);

% Calculate the phase
kf = 50;  % Frequency factor
phi_t = 2*pi*kf * cumsum(m_t) * Ts;  % Cumulative sum to get the phase

% Plot the time-domain representation of the phase
figure;
subplot(2, 1, 1);
plot(t, phi_t);
title('Phase phi(t)');
xlabel('t');
ylabel('phi(t)');
grid on;

% Calculate and plot the magnitude spectrum of the phase
n = length(phi_t);            % For normalize 
fshift = (-n/2:n/2-1)*(fs/n); % Calculate f parameter
Phi_f = abs(fftshift(fft(phi_t)))/length(phi_t);  % Calculate Phi(f)

subplot(2, 1, 2);
plot(fshift, Phi_f);
title('Magnitude Spectrum of Phase: Phi(f)');
xlabel('Frequency (Hz)');
ylabel('|Phi(f)|');
grid on;
xlim([-4, 4]);

% Define the modulated signal y(t)
y_t = 5*cos(500*pi*t + phi_t);   % y(t) = 5cos(2pi250t+ phi(t))

% Plot y(t) in the time domain
figure;
subplot(2, 1, 1);
plot(t, y_t);
title('Modulated Signal y(t)');
xlabel('t');
ylabel('y(t)');
grid on;
xlim([1.95, 2.05]);

% Calculate and plot the magnitude spectrum of y(t)
n_y = length(y_t);                    % For normalize 
fshift_y = (-n_y/2:n_y/2-1)*(fs/n_y); % Calculate f parameter
Y_f = abs(fftshift(fft(y_t)))/length(y_t); % Calculate Y(f)

subplot(2, 1, 2);
plot(fshift_y, Y_f);
title('Magnitude Spectrum of y(t): Y(f)');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
grid on;

% Differentiation
dy_dt = diff(y_t) / Ts;

% Plot the derivative in the time domain
figure;
subplot(2, 1, 1);
plot(t(2:end), dy_dt);
title('Derivative of y(t) in Time Domain');
xlabel('t');
ylabel('dy/dt');
grid on;
xlim([1.95, 2.05]);

% Calculate and plot the magnitude spectrum of dy(t)/dt
n_dy = length(dy_dt);                     % For normalize 
fshift_dy = (-n_dy/2:n_dy/2-1)*(fs/n_dy); % Calculate f parameter
dY_f = abs(fftshift(fft(dy_dt)))/length(dy_dt); 

subplot(2, 1, 2);
plot(fshift_dy, dY_f);
title('Magnitude Spectrum of dy(t)/dt: dY(f)');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
grid on;

% Rectification (Envelope Detection)
% Capacitor parameters
tau_d = 1e-2;   % Discharging time constant

% Initialize capacitor voltage
V_c = zeros(size(t));

V_c(1) = dy_dt(1); % Initial condition: Voltage at the first time point is the derivative of the signal
V_max = dy_dt(1);  % Initial maximum voltage
t0 = t(1);         % Initial time

for i = 2:length(dy_dt)
    % Check if the derivative is positive and greater than the previous voltage
    if dy_dt(i) > 0 && dy_dt(i)>V_c(i-1)
        % Charging phase
        if dy_dt(i) >= dy_dt(i-1)
            V_max = dy_dt(i);  % Update the maximum voltage during the charging phase
            t0 = t(i);         % Update the time when the maximum voltage occurs
            V_c(i) = dy_dt(i); % Set the capacitor voltage to the current derivative value
        % Discharging phase
        else
            % Exponential decrease during the discharging phase
            V_c(i) = V_max .* exp(-((t(i)-t0)/tau_d));
        end
    % Discharging phase
    else
        % Exponential decrease during the discharging phase
        V_c(i) = V_max .* exp(-((t(i)-t0)/tau_d));
    end
end

% Plot the envelope in the time domain
figure;
subplot(2, 1, 1);
plot(t(1:end), V_c);       % Adjusted the plotting indices
title('Envelope of y(t)');
xlabel('t');
ylabel('Envelope(t)');
grid on;
xlim([-2, 6]);

% Frequency domain representation of the envelope
n_envelope = length(V_c);
fshift_envelope = (-n_envelope/2:n_envelope/2-1)*(fs/n_envelope);
Envelope_f = fftshift(fft(V_c));

% Plot the magnitude spectrum of the envelope
subplot(2, 1, 2);
plot(fshift_envelope, abs(Envelope_f));
title('Magnitude Spectrum of Envelope: Y(f)');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
grid on;
xlim([-50, 50]);

% Apply a low-pass filter in the frequency domain
cutoff_frequency = 20;  % Cutoff_frequency
low_pass_filter = exp(-(fshift_envelope/cutoff_frequency).^2); % Gaussian LPF formula
low_pass_filter = low_pass_filter / max(low_pass_filter);  % Normalize the filter
Filtered_Envelope_f = (Envelope_f .* low_pass_filter);

% Plot the magnitude spectrum of the filtered envelope
figure;
subplot(2, 1, 1);
plot(fshift_envelope, abs(Filtered_Envelope_f));
title('Magnitude Spectrum of Filtered Envelope');
xlabel('Frequency (Hz)');
ylabel('|Filtered Envelope(f)|');
grid on;
xlim([-50, 50]);

% Inverse Fourier transform to obtain the filtered envelope in the time domain
filtered_envelope_t = ifft(ifftshift(Filtered_Envelope_f), 'symmetric');

% Plot the filtered envelope in the time domain
subplot(2, 1, 2);
plot(t(1:end), filtered_envelope_t);
title('Filtered Envelope in Time Domain');
xlabel('t');
ylabel('Filtered Envelope(t)');
grid on;

% DC filtering
dc_filtered_envelope_t = filtered_envelope_t - mean(filtered_envelope_t);

% Plot the DC-filtered envelope in the time domain
figure;
subplot(2, 1, 1);
plot(t(1:end), dc_filtered_envelope_t);
title('DC-Filtered Envelope in Time Domain');
xlabel('t');
ylabel('DC-Filtered Envelope(t)');
grid on;

% Frequency domain representation of the DC-filtered envelope
n_dc_filtered = length(dc_filtered_envelope_t);
fshift_dc_filtered = (-n_dc_filtered/2:n_dc_filtered/2-1)*(fs/n_dc_filtered);
DC_Filtered_Envelope_f = fftshift(fft(dc_filtered_envelope_t));

% Plot the magnitude spectrum of the DC-filtered envelope
subplot(2, 1, 2);
plot(fshift_dc_filtered, abs(DC_Filtered_Envelope_f));
title('Magnitude Spectrum of DC-Filtered Envelope');
xlabel('Frequency (Hz)');
ylabel('|DC-Filtered Envelope(f)|');
grid on;
xlim([-50, 50]);

% We multiply 2pikf and 5 (from carrier) so we must divide 5*2pikf
demodmt = dc_filtered_envelope_t / (5*2*pi*kf);
figure;
subplot(2, 1, 1);
plot(t(1:end), demodmt);
title('Demodulation Signal in Time Domain');
xlabel('t');
ylabel('Demod_m(t)');
grid on;

% Frequency domain representation of the demodulated signal
n_demod = length(demodmt);
fshift_demod = (-n_demod/2:n_demod/2-1)*(fs/n_demod);
Demodmt_f = fftshift(fft(demodmt)/fs);

% Plot the magnitude spectrum of the demodulated signal
subplot(2, 1, 2);
plot(fshift_demod, abs(Demodmt_f));
title('Magnitude Spectrum of Demodulation Signal');
xlabel('Frequency (Hz)');
ylabel('|Demod_m(f)|');
grid on;
xlim([-50, 50]);

% Plot m(t) and demod_m(t) together
figure;
subplot(2, 1, 1);
plot(t, m_t, 'b');  % Blue line for m(t)
hold on;
plot(t, demodmt, 'r');  % Red line for demod_m(t)
title('Signal m(t) and Demodulated Signal in Time Domain');
xlabel('t');
ylabel('Amplitude');
legend('m(t)', 'Demod_m(t)');
grid on;

% Plot the magnitude spectrum of m(t) and demod_m(t)
subplot(2, 1, 2);
plot(fshift, M_f, 'b');  % Blue line for M(f)
hold on;
plot(fshift_demod, abs(Demodmt_f), 'r');  % Red line for Demod_m(f)
title('Magnitude Spectrum of m(t) and Demodulated Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('|M(f)|', '|Demod_m(f)|');
grid on;
xlim([-50, 50]);