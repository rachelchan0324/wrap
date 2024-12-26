%% Simulation Setup

Fc = 1e6; % Carrier frequency in Hz
Fs = 4e6; % Sampling rate
Rs = 50e3; % Symbol rate in symbols/sec (baud)
B = 350e3; % Bandwidth

sps = Fs / Rs; % Number of samples per symbol
num_symbols = 296; % Number of symbols/bits to send
num_samples = num_symbols * sps; % Number of discrete samples

% Use these guys for plotting
t = linspace(0, num_symbols / Rs, num_samples); % Time variable
f = linspace(-Fs/2,Fs/2, num_samples); % Frequency variable

% Clear all previous figures before starting
close all;

%% Mapping Bits to Symbols

% Below are the bits for the binary message that we want to send

bits = [0 1 0 1 0 1 1 1 0 1 1 0 1 1 1 1 0 1 1 1 0 1 1 1 0 0 1 0 0 0 0 0 0 1,... 
        0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 1 0 0 1 1 0 0 0 1 1 0 1 1 1 1 0 1 1 1,... 
        0 1 1 0 0 1 1 0 0 1 0 1 0 0 1 0 0 0 0 0 0 1 0 1 0 1 1 1 0 1 0 1 0 0,... 
        1 0 0 1 0 0 0 0 0 1 0 1 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 1 0 0 1 1,... 
        0 1 1 0 1 1 1 1 0 0 1 0 0 0 0 0 0 1 0 0 1 1 0 1 0 1 1 1 0 1 0 1 0 1,... 
        1 0 0 0 1 1 0 1 1 0 1 0 0 0 0 0 1 0 0 0 0 0 0 1 0 1 0 1 0 0 0 1 1 0,... 
        1 0 0 0 0 1 1 0 1 0 0 1 0 1 1 1 0 0 1 1 0 0 1 0 0 0 0 0 0 1 0 0 1 0,... 
        0 1 0 1 1 1 0 0 1 1 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1 1 0 1 1 0 1 1 1 1,... 
        0 1 1 0 1 1 1 1 0 1 1 0 1 1 0 0 0 0 1 0 0 0 0 1];

% Convert to BPSK Symbols (1's and -1's)

% Initialize symbols array with same size as bits
symbols = zeros(size(bits));

% Loop through each bit and convert
for i = 1:length(bits)
    if bits(i) == 0
        symbols(i) = -1;
    else  % bits(i) == 1
        symbols(i) = 1;
    end
end

% Visualize BPSK Symbols in a Constellation Diagram
scatterplot(symbols);
title("Transmitted Symbols");

%% Represent Symbols as Time Shifted Deltas

% TODO 1.2.1: Use upsample to form a delta train
deltas = upsample(symbols, sps);

% Visualize Deltas (xlimited to only first ~20 symbols)
figure;
stem(t, deltas);
title("Upsampled Deltas");
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 0.0004]);

%% Filter Deltas with Transmitter Filter

% TODO 1.3.1: Generate rectangular pulse shape filter
% Hint: use the ones function
ps_filter = ones(1, sps);

% TODO 2.2: Replace the rect with an SRRC as the pulse shape filter 
% (ignore this if you are not at assignment section 2 yet)

% TODO 1.3.2: Convolve the deltas with the rectangular window
% Hint: use the conv function
transmitted_baseband = conv(deltas, ps_filter, 'same');

% Visualize Transmitted Baseband Signal
figure;
plot(t, transmitted_baseband);
title("Transmited Baseband Signal");
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 0.0004]);

%% Modulate Baseband Signal to Passband

% TODO 1.4.1: Modulate the baseband signal to passband
complex_carrier = exp(-1j*2*pi*Fc*t);  % e^(-j2πfct)
analytic_signal = transmitted_baseband .* complex_carrier;
transmitted_signal = real(analytic_signal);

fft_tx = fft(transmitted_signal);
% Shift and take magnitude
fft_mag = fftshift(abs(fft_tx));

% Now plot

% Create a figure with two subplots
figure;

% First subplot: FFT of transmitted signal
subplot(2,1,1);
plot(f, fft_mag);
title("FFT of Transmitted Signal");
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Second subplot: Time domain transmitted signal
subplot(2,1,2);
plot(t, transmitted_signal);
title('Transmitted Passband Signal');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 0.0004]); % Limit to see first few symbols clearly

%% Transmit Through Wireless Channel

% TODO 2.1: Apply the channel bandlimit
% (ignore this if you are not at assignment section 2 yet)

% TODO 1.5.1: Add some AWGN to the signal

received_signal = awgn(transmitted_signal, 20, 'measured');

% TODO 1.5.2: Plot the Corrupted Signal

% Plot frequency domain of the received (noisy) signal
figure;
% Calculate FFT and shift to center
fft_rx = fft(received_signal);
fft_mag = fftshift(abs(fft_rx));

% Make sure frequency vector matches length of FFT
f = linspace(-Fs/2, Fs/2, length(received_signal));

% Plot
subplot(2,1,1);  % Time domain plot
plot(t, received_signal);
title('Received Signal with AWGN (Time Domain)');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 0.0004]);

subplot(2,1,2);  % Frequency domain plot
plot(f, fft_mag);
title('Received Signal with AWGN (Frequency Domain)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%% Demodulation

% TODO 1.6.1: Start demodulation by mixing the received signal by mixing with the sinusoids
I_ = received_signal .* cos(2*pi*Fc*t);  % Mix with cosine (not 2*cos)
Q_ = -received_signal .* sin(2*pi*Fc*t); % Mix with negative sine (not -2*sin)

% Lowpass filter the previous results to get I and Q
I = lowpass(I_, 1.5e6, Fs);
Q = lowpass(Q_, 1.5e6, Fs);

% TODO 1.6.2: Plot the In Phase Component
figure;
plot(t, I);  % Plot the filtered I component versus time
title('Demodulated BPSK In-Phase Component');
xlabel('Time');
ylabel('Amplitude');
grid on;

%% Apply Receiver Impulse Response

% TODO 1.7.1: Combine I and Q components to construct the received baseband signal (I + jQ)
received_baseband = I + 1j * Q

% TODO 1.7.2: Convolve the received baseband signal with the receiver filter
% received_samples = conv(received_baseband, , 'same');

%% Sample and Detect Symbols

% TODO 1.8.1: (Naively) sample the received signal at the symbol rate to get the received symbols
% Hint: Use the downsample function
received_symbols = downsample(received_baseband, sps, 1);  % Sample every sps samples, starting at 1

% TODO 1.8.2: Visualize the received symbols in a constellation diagram (scatterplot)
figure;
scatter(real(received_symbols), imag(received_symbols), 'b.');
title('BPSK Constellation Diagram');
xlabel('In-Phase (I)');
ylabel('Quadrature (Q)');
grid on;
axis equal;  % Make the plot square
hold on;
% Add reference points for ideal constellation
plot([-1 1], [0 0], 'rx', 'MarkerSize', 10);  % Show ideal constellation points

% TODO 1.8.3: Threshold detection
received_symbols = sign(real(received_symbols));  % Convert to ±1 using sign function

%% Map Symbols back to Bits

% TODO 1.9.1: Map the detected symbols back to bits (turn them back to 1's
% and 0's
detected_bits = (received_symbols + 1) / 2;  % Converts -1 → 0 and +1 → 1

% Calculate Bit Error Rate
num_bit_errors = sum(detected_bits(1:length(bits)) ~= bits)
BER = 100 * num_bit_errors / length(bits)

%% Decode and Display Message

% If you have successfully decoded the message, you should see the message "Congrats! Simulation is complete"
s = string(detected_bits);
s = strjoin(s);
s = strrep(s, " ", "");
disp(s);
s_len = length(char(s))
inputString = char(s);
binaryString = inputString(1:end-mod(length(inputString),8));
binaryChunks = reshape(binaryString, 8, []).';
asciiChars = char(bin2dec(binaryChunks)).';
disp(asciiChars);