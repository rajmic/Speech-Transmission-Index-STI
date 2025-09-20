function plotIRVisualization(audiodata, titleStr)
% PLOTIRVISUALIZATION visualizes exponential sweep and inverse filter waveforms
% and spectrograms for IR-based STI measurement.
%
%   plotIRVisualization(AUDIODATA) uses AUDIODATA struct with fields:
%     .audio  — exponential sweep signal,
%     .audio2 — inverse filter,
%     .fs     — sampling rate in Hz.
%
%   plotIRVisualization(AUDIODATA, TITLESTR) prefixes subplot titles with TITLESTR.
%
%   Generates two figures:
%     Figure 1: time-domain waveforms (sweep & inverse filter)
%     Figure 2: spectrograms (sweep & inverse filter)
%
% Copyright Šimon Cieslar, Brno University of Technology, 2024-2025

    if nargin < 2
        titleStr = '';
    end
    
    % Extract signals and sampling rate
    sweep = audiodata.audio;
    invF  = audiodata.audio2;
    fs    = audiodata.fs;

    % Time vectors
    t_sweep = (0:length(sweep)-1) / fs;
    t_inv   = (0:length(invF)-1)  / fs;

    %--- Figure 1: Time-domain waveforms ---
    figure;
    subplot(2,1,1);
    plot(t_sweep, sweep);
    xlabel('Time (s)'); ylabel('Amplitude');
    if isempty(titleStr)
        title('Exponential Sweep');
    else
        title([titleStr ' - Exponential Sweep']);
    end
    grid on;

    subplot(2,1,2);
    plot(t_inv, invF);
    xlabel('Time (s)'); ylabel('Amplitude');
    if isempty(titleStr)
        title('Inverse Filter');
    else
        title([titleStr ' - Inverse Filter']);
    end
    grid on;

    %--- Figure 2: Spectrograms ---
    figure;
    subplot(2,1,1);
    spectrogram(sweep, hamming(1024), 512, 1024, fs, 'yaxis');
    if isempty(titleStr)
        title('Spectrogram: Exponential Sweep');
    else
        title([titleStr ' - Spectrogram: Exponential Sweep']);
    end
    ylabel('Frequency (kHz)'); xlabel('Time (s)'); colorbar;

    subplot(2,1,2);
    spectrogram(invF, hamming(1024), 512, 1024, fs, 'yaxis');
    if isempty(titleStr)
        title('Spectrogram: Inverse Filter');
    else
        title([titleStr ' - Spectrogram: Inverse Filter']);
    end
    ylabel('Frequency (kHz)'); xlabel('Time (s)'); colorbar;
end