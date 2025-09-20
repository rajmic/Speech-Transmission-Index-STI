%% Demonstration file for IR STI Evaluation
% This script serves as a demonstration for generating a test signal, applying a mono reverb effect,
% and visualizing both the original and processed signals for IR STI evaluation.
clear; close all; clc;

%% Generate Exponential Sweep Test Signal
% Define parameters for the exponential sweep test signal.
dur = 2;         % Duration of the test signal in seconds
start_freq = 20;  % Starting frequency in Hz
end_freq = 20000; % Ending frequency in Hz
fs = 48000;       % Sampling frequency in Hz
reverse = 0;      % Do not reverse the signal (0 = original order)
rcos_ms = 15;     % Fade-in and fade-out duration in milliseconds

% Generate the test signal (exponential sweep and inverse signal)
audiodata = IR_signal_exp_sweep(dur, start_freq, end_freq, fs, reverse, rcos_ms);
fprintf('IR STI test signal has been generated and saved in IR_signal.mat.\n');

%% Create Time Vectors for Signal Visualization
% Construct time vectors for plotting the original and processed signals in the time domain.
t_audio   = (0:length(audiodata.audio)-1) / audiodata.fs;
t_audio2  = (0:length(audiodata.audio2)-1) / audiodata.fs;

%% Plot Original Test Signals
% Plot the exponential sweep and its inverse in separate subplots for visual analysis.
figure;

% Plot exponential sweep (audio1)
subplot(3,1,1);
plot(t_audio, audiodata.audio);
xlabel('Time (s)');
ylabel('Amplitude');
title('Exponential Sweep (audio1)');
grid on;

% Plot inverse signal (audio2)
subplot(3,1,2);
plot(t_audio2, audiodata.audio2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Inverse Signal (audio2)');
grid on;

%% Apply REVERB Effect Using Audio Toolbox Reverberator
% Create a reverberator object with specified parameters to simulate a room reverb effect.
reverb = reverberator('PreDelay', 0.01, 'SampleRate', audiodata.fs, 'WetDryMix', 0.8, 'DecayFactor', 0.2);

% Convert the exponential sweep signal to mono if necessary
if size(audiodata.audio,2) > 1
    monoSignal = mean(audiodata.audio, 2);  % Average channels for a stereo signal
else
    monoSignal = audiodata.audio;
end

% Apply the reverb effect to the mono signal
reverbSignal = reverb(monoSignal);

% Ensure the processed signal is mono by averaging channels if needed
if size(reverbSignal,2) > 1
    reverbSignal = mean(reverbSignal, 2);
end

% Create a time vector for the reverberated signal
t_reverb = (0:length(reverbSignal)-1) / audiodata.fs;

% Save the processed (reverberated) signal to a WAV file
OutputFileName = 'IR_Signal_exp_sweep_with_reverb.wav';
audiowrite(OutputFileName, reverbSignal, audiodata.fs, 'BitsPerSample', 24);
fprintf('Reverberated mono signal has been saved as %s.\n', OutputFileName);

% Plot the reverberated signal for visual inspection
subplot(3,1,3);
plot(t_reverb, reverbSignal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Exponential Sweep with Reverb (audio1) - Mono');
grid on;

%% Convolution Using Built-in conv Function
% Perform convolution of the reverberated signal with the inverse signal using MATLAB's built-in conv function.
conv_sig = conv(reverbSignal, audiodata.audio2);

% Create a time vector for the convolved signal
t_conv = (0:length(conv_sig)-1) / audiodata.fs;

% Plot the convolved signal to observe the resulting impulse response segment
figure;
plot(t_conv, conv_sig);
xlabel('Time (s)');
title('IR - Convolution of Signals (Built-in conv Function)');
grid on;

%% Convolution of CLEAN GENERATED SIGNALS Using Built-in conv Function
% Perform convolution of the reverberated signal with the inverse signal using MATLAB's built-in conv function.
conv_sig_clean = conv(audiodata.audio, audiodata.audio2);

% Create a time vector for the convolved signal
t_conv_clean = (0:length(conv_sig_clean)-1) / audiodata.fs;

% Plot the convolved signal to observe the resulting impulse response segment
figure;
plot(t_conv_clean, conv_sig_clean);
xlabel('Time (s)');
title('IR - Convolution of CLEAN GENERATED Signals (Built-in conv Function)');
grid on;

%% Convolution Using Custom Implementation
% Alternatively, perform convolution using a custom convolution function to verify the results.
conv_signal = convolution(reverbSignal, audiodata.audio2, fs);

% Create a time vector for the custom convolved signal
t_conv2 = (0:length(conv_signal)-1) / audiodata.fs;

% Plot the convolved signal from the custom implementation for comparison
figure;
plot(t_conv2, conv_signal);
xlabel('Time (s)');
title('IR - Convolution of Signals (Custom Implementation)');
grid on;

%% Plot Combined Frequency Responses
% Compute FFT parameters for the exponential sweep, inverse filter, and impulse response.

% Exponential Sweep
NFFT_sweep = 2^nextpow2(length(audiodata.audio));
f_sweep = linspace(0, audiodata.fs/2, NFFT_sweep/2+1);
S_fft = fft(audiodata.audio, NFFT_sweep);
S_mag = 20*log10(abs(S_fft(1:NFFT_sweep/2+1)));

% Inverse Filter
NFFT_inv = 2^nextpow2(length(audiodata.audio2));
f_inv = linspace(0, audiodata.fs/2, NFFT_inv/2+1);
Sinv_fft = fft(audiodata.audio2, NFFT_inv);
Sinv_mag = 20*log10(abs(Sinv_fft(1:NFFT_inv/2+1)));

% Impulse Response (from conv_sig)
NFFT_ir = 2^nextpow2(length(conv_sig_clean));
f_ir = linspace(0, audiodata.fs/2, NFFT_ir/2+1);
IR_fft = fft(conv_sig_clean, NFFT_ir);
IR_mag = 20*log10(abs(IR_fft(1:NFFT_ir/2+1)));

% Plot all magnitude responses in a single graph
figure;
semilogx(f_sweep, S_mag, 'b', 'LineWidth', 1.5); hold on;
semilogx(f_inv, Sinv_mag, 'r', 'LineWidth', 1.5);
semilogx(f_ir, IR_mag, 'g', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Combined Frequency Responses');
legend('Exponential Sweep', 'Inverse Filter', 'Impulse Response');
grid on;
hold off;

%% Plot Spectrograms of Test Signals
% Plotting the spectrograms for the exponential sweep and its inverse signal.
figure;
subplot(2,1,1);
spectrogram(audiodata.audio, hamming(1024), 512, 1024, audiodata.fs, 'yaxis');
title('Spectrogram of Exponential Sweep (audio1)');
xlabel('Time (s)');
ylabel('Frequency (kHz)');
colorbar;

subplot(2,1,2);
spectrogram(audiodata.audio2, hamming(1024), 512, 1024, audiodata.fs, 'yaxis');
title('Spectrogram of Inverse Signal (audio2)');
xlabel('Time (s)');
ylabel('Frequency (kHz)');
colorbar;

%% Calculation MTF and STI
[STI, mk, STIPA_IR] = STI_IR(conv_signal, fs, 'doTable', 1);
fprintf('Computed STI value of signal: %.2f.\n', STI)
fprintf('Computed STI - STIPA(IR) value of signal: %.2f.\n', STIPA_IR)

%% Compute STI value adjusted to auditory masking and ambient noise
Lsk = [72.2, 72.3, 70.1, 59.7, 51.5, 42.8, 36.4]; % measured signal levels
Lnk = [41.7, 42.0, 44.3, 38.1, 24.2, 21.0, 19.6]; % measured levels of ambient noise
STI_ = STI_IR(conv_signal, fs, 'Lsk', Lsk, 'Lnk', Lnk);
fprintf(['Computed STI value adjusted for auditory masking ' ...
    'and ambient noise: %.2f.\n'], STI_)

%% Plot Test Signal, Inverse Filter, Waveform, and Spectrogram
figure('Position',[100, 100, 1000, 800]);
% Top-left: Exponential Sine-Sweep Test Signal
subplot(2,1,1);
plot(t_audio, audiodata.audio, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude'); xlim([0 1.2]);
title('Waveform of Exponential Sine-Sweep Test Signal');
grid on;

% Bottom-left: Inverse Filter for Exponential Sweep
subplot(2,1,2);
spectrogram(audiodata.audio, 256, 250, 256, audiodata.fs, 'yaxis');
title('Spectrogram of Generated Test Signal');
xlabel('Time (s)');
ylabel('Frequency (kHz)');
grid on;

figure('Position',[100, 100, 1000, 800]);
% Top-right: Waveform of Generated Test Signal
subplot(2,1,1);
plot(t_audio2, audiodata.audio2, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude'); xlim([0 1.2]);
title('Waveform of Inverse Filter for Exponential Sine-Sweep');
grid on;

% Bottom-right: Spectrogram of Generated Test Signal
subplot(2,1,2);
spectrogram(audiodata.audio2, 256, 250, 256, audiodata.fs, 'yaxis');
title('Spectrogram of Generated Inverse Filter');
xlabel('Time (s)');
ylabel('Frequency (kHz)');
grid on;