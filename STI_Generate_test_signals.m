%% Generate all STI test signals
% This MATLAB script generates a variety of standardized test signals for
% Speech Transmission Index (STI) measurements.
%
% It produces:
%   1. Full-band STI signals with configurable segment durations and silence gaps
%   2. STIPA signals of varying lengths for rapid STI assessment
%   3. Exponential sweep signals (ESS) and inverse filters for impulse-response (IR) measurement
%
% All signals are sampled at 48 kHz and saved as 24-bit WAV files into a
% dedicated output folder.
clear; close all; clc;

%% PARAMETERS
fs                   = 48000;             % Sampling frequency (Hz)
outputFolder         = 'test_signals';    % Folder to store generated files (.wav, .mat)
fullDurations        = [12, 12];          % Full‑STI segment lengths (s)
fullSilences         = [0, 1];            % Silence gaps between Full‑STI segments (s)
stipaDurations       = 25;                % STIPA signal lengths (s)
irDurations          = [2, 5];            % ESS sweep lengths (s)
irFreqRanges         = {[20, 20000], ...  % Sweep start/end frequencies (Hz)
                        [20, 10000]};

%% SETUP OUTPUT FOLDER
% Create output directory for all generated signals
if ~exist(outputFolder, 'dir')             % Check if the folder already exists
    mkdir(outputFolder);                   % Create the folder if it doesn't exist
end

%% FULL STI Signals Generation
% Generates concatenated narrowband noise segments (98 bands) with
% configurable duration and inter‑segment silence for full‑band STI.

for idx = 1:numel(fullDurations)
    dur    = fullDurations(idx);
    sil    = fullSilences(idx);
    sig    = generateFullSTISignal(dur, fs, sil, 'doPlot', 1);
    fprintf(['Generated %d segments of %d seconds with silence gap of %d ' ...
        'seconds between each segments, sampled at %d Hz.\n'], 98, dur, sil, fs);

    fname  = sprintf('Full_STI_Signal_%ds_%ds.wav', dur, sil);
    audiowrite(fullfile(outputFolder, fname), sig, fs, 'BitsPerSample', 24);
end

%% STIPA Signals Generation
% STIPA (Speech Transmission Index for Public Address systems) uses a
% single rapidly modulated signal covering multiple octave bands for
% quick STI estimation in live venues.

for dur = stipaDurations
    sig = generateStipaSignal(dur, fs, 'doPlot', 1);
    fprintf('Generated %g seconds of STIPA test signal sampled at %d Hz.\n', dur, fs);

    fname = sprintf('STIPA_Signal_%ds.wav', dur);
    audiowrite(fullfile(outputFolder, fname), sig, fs, 'BitsPerSample', 24);
end

%% IR STI - Exponential Sweep (ESS) & Inverse Filter Generation
% Exponential sine sweep signals are used to measure the impulse response
% of an acoustic system. The generated IR.mat file contains both the
% sweep and its inverse filter for deconvolution.

for dur = irDurations
    for r = 1:numel(irFreqRanges)
        fr    = irFreqRanges{r};
        f0    = fr(1);
        f1    = fr(2);
        ir    = IR_signal_exp_sweep(dur, f0, f1, fs, 'doPlot', 1);
        fprintf('Generated %g seconds of IR ESS sampled at %d Hz with frequency range %g - %d Hz.\n', dur, fs, f0, f1);

        % Filenames
        tag      = sprintf('%ds_%dHz-%dkHz', dur, f0, f1/1000);
        fname_s  = fullfile(outputFolder, ['IR_ESS_' tag '.wav']);
        fname_i  = fullfile(outputFolder, ['IR_ESS_' tag '_inverse.wav']);
        fname_m  = fullfile(outputFolder, ['IR_ESS_' tag '.mat']);

        % Write audio .wav and .mat
        audiowrite(fname_s, ir.audio,  ir.fs, 'BitsPerSample', 24);
        audiowrite(fname_i, ir.audio2, ir.fs, 'BitsPerSample', 24);
        save(      fname_m, 'ir');
    end
end

%% OUTPUT SUMMARY
%
% Files in “Test signals/”:
%   Full_STI_Signal_<dur>s_<silence>s.wav
%   STIPA_Signal_<dur>s.wav
%   IR_ESS_<dur>s_<start>Hz-<end>kHz.wav
%   IR_ESS_<dur>s_<start>Hz-<end>kHz_inverse.wav
%   IR_ESS_<dur>s_<start>Hz-<end>kHz.mat