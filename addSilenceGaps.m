function modifiedSignal = addSilenceGaps(signal, fs, silenceStartDuration, silenceEndDuration)
% ADDSILENCEGAPS Adds silence at the beginning and end of an audio signal.
%
%   MODIFIEDSIGNAL = ADDSILENCEGAPS(SIGNAL, FS, SILENCESTARTDURATION, SILENCEENDDURATION)
%   inserts SILENCESTARTDURATION seconds of silence before SIGNAL and
%   SILENCEENDDURATION seconds of silence after SIGNAL. FS is in Hz.
%
%   Inputs:
%       signal              - Audio signal (vector or matrix)
%       fs                  - Sampling frequency (in Hz)
%       silenceStartDuration- Duration of silence to add at the beginning (in seconds)
%       silenceEndDuration  - Duration of silence to add at the end (in seconds)
%
%   Output:
%       modifiedSignal      - The audio signal with silence added at the beginning and end.
%
%   Example:
%       y_mod = addSilenceGaps(y, fs, 2, 3);
%
% Copyright Šimon Cieslar, Brno University of Technology, 2024-2025

% Calculate the number of samples for the silence durations
silenceStartSamples = round(silenceStartDuration * fs);
silenceEndSamples   = round(silenceEndDuration * fs);

% Create silence segments (zero vectors or matrices if multi-channel)
silenceStart = zeros(silenceStartSamples, size(signal, 2));
silenceEnd   = zeros(silenceEndSamples, size(signal, 2));

% Concatenate the silence at the beginning and at the end of the signal
modifiedSignal = [silenceStart; signal; silenceEnd];

end