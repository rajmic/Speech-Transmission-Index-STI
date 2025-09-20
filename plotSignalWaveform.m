function plotSignalWaveform(signal, fs, titleStr)
% PLOTSIGNALWAVEFORM plots a time-domain waveform of any audio signal.
%
%   plotSignalWaveform(SIGNAL, FS) displays SIGNAL over 0:(N−1)/FS seconds,
%   labels axes, adds grid and default title “Signal Waveform.”
%
%   plotSignalWaveform(SIGNAL, FS, TITLESTR) uses TITLESTR as the plot title.
%
%   Inputs:
%     SIGNAL   — [N×1] audio sample vector.
%     FS       — sampling frequency in Hz.
%     TITLESTR — (optional) plot title (default: 'Signal Waveform').
%
%   This helper is used to visualize time-domain signals generated for
%   Full STI and STIPA measurements.
%
% Copyright Šimon Cieslar, Brno University of Technology, 2024-2025

    if nargin < 3
        titleStr = 'Signal Waveform';
    end
    t = (0:length(signal)-1) / fs;
    figure;
    plot(t, signal);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(titleStr);
    grid on;
    
end