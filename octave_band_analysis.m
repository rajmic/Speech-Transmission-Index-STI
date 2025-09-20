% OCTAVE_BAND_ANALYSIS computes octave-band and A-weighted levels from a recording.
%
%   This script executes the following steps:
%     1) Calibration using a calibrator recording (94 dB SPL = 1 Pa)
%     2) Octave-band filtering via octaveFilter
%     3) RMS calculation and conversion to dB SPL
%     4) Apply A-weighting offsets for each band
%
%   Requirements:
%     • Signal Processing Toolbox (octaveFilter)
%     • Audio Toolbox (audioread)
%
%   Output:
%     Table with variables:
%       CenterFreq_Hz    Level_dB_SPL    Level_dB_A_SPL
%
% Copyright Šimon Cieslar, Brno University of Technology, 2024-2025.

% Load calibrator recording
[calFile, calPath] = uigetfile('*.wav', 'Select calibrator file (94 dB SPL)');
if isequal(calFile, 0)
    error('No calibrator file selected.');
end
[calData, fs] = audioread(fullfile(calPath, calFile));
if size(calData,2) > 1
    calData = mean(calData, 2);
end

% Compute calibration constant (94 dB SPL = 1 Pa)
rmsCal = rms(calData);
pCal   = 1;  % Pa
K      = pCal / rmsCal;

% Load measurement recording
[measFile, measPath] = uigetfile('*.wav', 'Select measurement file');
if isequal(measFile, 0)
    error('No measurement file selected.');
end
[measData, fs2] = audioread(fullfile(measPath, measFile));
if fs2 ~= fs
    error('Sampling rate mismatch between calibrator and measurement.');
end
if size(measData,2) > 1
    measData = mean(measData, 2);
end

% Define octave-band center frequencies and A-weighting offsets
fc         = [125, 250, 500, 1000, 2000, 4000, 8000];
Aweights   = [-16.1, -8.6, -3.2, 0.0, 1.2, 1.0, -1.1];
filterOrder = 8;
nBands      = numel(fc);

% Preallocate results
L_oct   = nan(nBands,1);
L_Aoct  = nan(nBands,1);

% Filter into octave bands and compute both levels
for k = 1:nBands
    % 1-octave filter for center frequency
    octFilt = octaveFilter(fc(k), '1 octave', 'SampleRate', fs, ...
                           'FilterOrder', filterOrder);

    % Filter raw digital signal
    bandSig = octFilt(measData);

    % Compute unweighted level (dB SPL)
    prms      = K * rms(bandSig);
    L_oct(k) = 20*log10(prms / 20e-6);

    % Compute A-weighted level (dB(A)) by applying offsets
    L_Aoct(k)= L_oct(k) + Aweights(k);
end

% Display results in a combined table
T = table(fc(:), L_oct, L_Aoct, ...
    'VariableNames', {'CenterFreq_Hz','Level_dB_SPL','Level_dB_A_SPL'});
fprintf('Octave-band levels (dB SPL) and A-weighted levels (dB(A)) for centers %s Hz:\n', mat2str(fc));
disp(T);