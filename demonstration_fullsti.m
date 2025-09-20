%% Demonstration file for Full STI Evaluation
% This script serves as a simple demonstration of FULL STI implementation.
clear; close all; clc;

%% Generate FULL STI test signal
duration = 10;       % Duration of individual signals in seconds
silenceDuration = 1; % Gap between individual signals in seconds
fs = 48000;          % Sample rate in Hz

FullSTISignal = generateFullSTISignal(duration, fs, silenceDuration, 'doPlot', 1);
fprintf(['A Full STI test signal consisting of 98 individual %g-second signals has been generated, ' ...
    'sampled at %d Hz, with a silence gap of %g seconds between each segments.\n'], duration, fs, silenceDuration)

% Save the final signal to the root directory
filename = 'Full_STI_Signal.wav';
audiowrite(filename, FullSTISignal, fs, 'BitsPerSample', 24);

%% Add silence gaps using the addSilenceGaps function
% This section adds a silence gap at the beginning and at the end of the signal.
silenceStartDuration = 2; % Duration of silence at the beginning (in seconds)
silenceEndDuration   = 0; % Duration of silence at the end (in seconds)

SignalWithSilenceGaps = addSilenceGaps(FullSTISignal, fs, silenceStartDuration, silenceEndDuration);

%% Apply the REVERB effect to the generated test signal, simulating its passage through a transmission channel, using the Audio Toolbox reverberator

% Create a reverberator object
reverb = reverberator('PreDelay', 0.01, 'SampleRate', fs, 'WetDryMix', 0.8, 'DecayFactor', 0.2); % Example settings

% Apply the reverb effect to the signal
reverbSignal = reverb(FullSTISignal);

% Convert stereo signal to mono if necessary
if size(reverbSignal, 2) > 1
    reverbSignal = mean(reverbSignal, 2); % Take the average of the two channels
end

% Save the processed signal
OutputFileName = 'Full_STI_Signal_with_reverb.wav';
audiowrite(OutputFileName, reverbSignal, fs, 'BitsPerSample', 24);

%% Compute STI value of generated test signal

STI = fullsti(FullSTISignal, fs, 'SegmentDuration', duration, 'SilenceDuration', silenceDuration, 'doPlot', 1, 'doTable', 1, 'SignalStart', 0);
fprintf('Computed STI value of generated test signal: %.2f.\n', STI)

%% Compute STI value of generated test signal with added silence gaps using the addSilenceGaps function
% Note: Automatic detection using an energy threshold is employed to determine the signal start.

STI = fullsti(SignalWithSilenceGaps, fs, 'SegmentDuration', duration, 'SilenceDuration', silenceDuration, 'doPlot', 1, 'doTable', 0);
fprintf('Computed STI value of generated test signal with added silence gaps: %.2f.\n', STI)

%% Compute STI value of generated test signal adjusted to auditory masking and ambient noise

Lsk = [72.2, 72.3, 70.1, 59.7, 51.5, 42.8, 36.4]; % Measured signal levels
Lnk = [41.7, 42.0, 44.3, 38.1, 24.2, 21.0, 19.6]; % Measured levels of ambient noise
STI = fullsti(FullSTISignal, fs, 'Lsk', Lsk, 'Lnk', Lnk, 'SegmentDuration', duration, 'SilenceDuration', silenceDuration, 'doTable', 0, 'doPlot', 1);
fprintf(['Computed STI value of generated test signal adjusted for auditory masking ' ...
  'and ambient noise: %.2f.\n'], STI)

%% Compute STI value of generated test signal with using the REFERENCE signal 
% and its sampling frequency FSREF. The STI value is assumed to be 1.

STI = fullsti(FullSTISignal, fs, FullSTISignal, fs, 'SegmentDuration', duration, 'SilenceDuration', silenceDuration, 'doPlot', 1);
fprintf('Computed STI value of generated test signal with using the REFERENCE signal and its sampling frequency FSREF.: %.2f.\n', STI)

%% Compute STI value with applied reverb and using the REFERENCE signal
 
STI_reverb = fullsti(reverbSignal, fs, FullSTISignal, fs,'SegmentDuration', duration, 'SilenceDuration', silenceDuration);
fprintf('Computed STI value with applied reverb and using the REFERENCE signal: %.2f.\n', STI_reverb)

%% Load recorded FULL STI signal after passing through the transmission channel

filenameRec = 'FullSti_10s_0s_silence_front_1.wav';
[FullStiRec, fsRec] = audioread(filenameRec);
fprintf('Loading measurement ''%s''.\n', filenameRec)
 
%% Compute STI value and using the REFERENCE signal

STI_rec = fullsti(FullStiRec, fsRec, FullSTISignal, fs);
fprintf('Computed STI value and using the REFERENCE signal: %.2f.\n', STI_rec)