function stipa_signal = generateStipaSignal(duration, varargin)
% GENERATESTIPASIGNAL generates the STIPA test signal for measurement of
% the Speech Transmission Index (STI) according to IEC 60268-16 (ed.3)
% standard using the direct STIPA (Speech Transmission Index for Public
% Addressing Systems) method.
%
%   STIPA_SIGNAL = GENERATESTIPASIGNAL(DURATION) generates DURATION seconds of
%   the STIPA test signal. Default sampling frequency is 96000 Hz, and the RMS
%   level is normalized to 0.07.
%
%   STIPA_SIGNAL = GENERATESTIPASIGNAL(DURATION, FS) generates DURATION seconds
%   of the STIPA test signal with a sampling frequency of FS Hertz.
%   The minimum FS is 22050 Hz.
%
%   STIPA_SIGNAL = GENERATESTIPASIGNAL(DURATION, FS, 'doPlot', DOPLOT)
%   when DOPLOT = 1 (default) displays the signal waveform via plotSignalWaveform,
%   and when DOPLOT = 0 suppresses plotting.
%
% References:
%   EN IEC 60268-16:2020 Sound system equipment — Part 16: Objective rating of speech
%   intelligibility by speech transmission index.
%   Záviška, P., Rajmic, P., Schimmel, J. (2024). Matlab Implementation 
%   of STIPA (Speech Transmission Index for Public Address Systems).
%
% Copyright Pavel Záviška, Brno University of Technology, 2023-2024

% Check number of input arguments
narginchk(1, 4);

    %% Parse input arguments
    [fs, doPlot] = parseInputs(duration, varargin{:});
    
    %% Generate pink noise
    N = duration * fs;
    pinkNoise = pinknoise(N);
    
    % Filtrate the pink noise into octave bands
    filteredPinkNoise = NaN(N, 7);
    filterOrder = 20;
    octaveBands = [125, 250, 500, 1000, 2000, 4000, 8000];
    
    for bandIdx = 1:length(octaveBands)
        octFilt = octaveFilter(octaveBands(bandIdx), '1/2 octave', ...
            'SampleRate', fs, 'FilterOrder', filterOrder);
        filteredPinkNoise(:, bandIdx) = octFilt(pinkNoise);
    end
    
    % Modulate the frequencies
    fm = [1.6, 1, 0.63, 2, 1.25, 0.8, 2.5; ... % modulation frequencies in Hz
        8, 5, 3.15, 10, 6.25, 4, 12.5];
    
    t = linspace(0, duration, N).';
    modulation = NaN(N, 7);
    
    for bandIdx = 1:length(octaveBands)
        modulation(:, bandIdx) = sqrt(0.5 * (1 + 0.55 * ...
            (sin(2 * pi * fm(1, bandIdx) * t) - ...
             sin(2 * pi * fm(2, bandIdx) * t))));
    end
    
    % Set levels of the octave bands
    levels = [-2.5, 0.5, 0, -6, -12, -18, -24]; % revision 5 band levels (dB)
    G = 10 .^ (levels / 20); % acoustic pressure of bands
    
    % Compute the final STIPA test signal
    stipa_signal = sum(filteredPinkNoise .* modulation .* G, 2);
    
    % Normalize the RMS of the final singal
    targetRMS = 0.07; % empirically derived from the character of the STIPA test signal
    stipa_signal = stipa_signal * targetRMS / rms(stipa_signal);

    %% Plot signal waveform if requested
    if doPlot
        titleStr = sprintf('STIPA: %ds at %d Hz', duration, fs);
        plotSignalWaveform(stipa_signal, fs, titleStr);
    end
end

%% Helper Functions

function [fs, doPlot] = parseInputs(duration, varargin)
% PARSEINPUTS Parses input arguments and returns sampling frequency,
% and plot option.
    p = inputParser;
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
    validSamplingFreq = @(x) isnumeric(x) && isscalar(x) && (x >= 22050);
    validPlotOption   = @(x) isnumeric(x) && isscalar(x) && (x == 0 || x == 1);
    
    addRequired(p, 'duration', validScalarPosNum);
    
    defaultFs = 96000;
    addOptional(p, 'fs', defaultFs, validSamplingFreq);
    
    defaultDoPlot = 1;
    addOptional(p, 'doPlot', defaultDoPlot, validPlotOption);

    parse(p, duration, varargin{:});
    fs     = p.Results.fs;
    doPlot = p.Results.doPlot;
end