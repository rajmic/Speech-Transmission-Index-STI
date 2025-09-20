function IR = convolution(chan_Signal, invTestSignal, fs, varargin)
% CONVOLUTION computes the channel impulse response by convolving the
% transmitted signal with an inverse test signal using FFT-based
% convolution, with optional trimming and scaling.
%
%   IR = CONVOLUTION(CHAN_SIGNAL, INVTESTSIGNAL, FS) returns the impulse
%   response IR obtained by linearly convolving CHAN_SIGNAL — the signal
%   that has passed through the transmission channel — with INVTESTSIGNAL —
%   the inverse of the test signal. Sampling frequency FS is in Hz.
%
%   IR = CONVOLUTION(..., 'SignalStart', T0) specifies the start time T0
%   (in seconds) for trimming leading silence/noise (default: automatic).
%
%   IR = CONVOLUTION(..., 'AutoTrim', TF) enables (TF=1) or disables
%   (TF=0) automatic trimming of leading silence (default: 1).
%
%   IR = CONVOLUTION(..., 'ScalingFactor', SCALEF) multiplies the
%   resulting IR by SCALEF (default: 1). If INVTESTSIGNAL is a struct with
%   field .IRscalingfactor, that value overrides SCALEF.
%
%   Note: if invTestSignal is a struct containing field
%     .IRscalingfactor, that value will be used as ScalingFactor
%     (overriding any other setting).
%
%   Inputs:
%     chan_Signal   – Column vector or matrix (each column a channel) of the
%                     measured output signal from the transmission channel.
%     invTestSignal – Column vector or matrix of the inverse test signal
%                     (deconvolution kernel), dimensioned to match chan_Signal.
%     fs            – Sampling frequency in Hz of chan_Signal and
%                     invTestSignal.
%     
%   Name-value Inputs:
%     'SignalStart'   – scalar start time in seconds (default: NaN = auto)
%     'AutoTrim'      – 1/0 enable automatic or manual trimming (default: 1)
%     'ScalingFactor' – scalar to multiply IR by (default: 1)
%
%   Output:
%     IR            – Column vector or matrix containing the estimated impulse
%                     response(s) of the channel. Length is length(chan_Signal)
%                     + length(invTestSignal) − 1.
%
% This function performs linear convolution using spectral multiplication.
% It ensures that the dimensions of the two input signals match and uses FFT-based
% convolution, which inherently applies zero-padding to achieve linear convolution.
%
%   Reference:
%       Alan V Oppenheim Ronald W. Schafer (2010). Discrete-Time Signal Processing, 3rd edition
%       [Available at: https://file.fouladi.ir/courses/dsp/books/%28Prentice-Hall%20Signal%20Processing%20Series%29%20Alan%20V.%20Oppenheim%2C%20Ronald%20W.%20Schafer-Discrete-Time%20Signal%20Processing-Prentice%20Hall%20%282009%29.pdf]
%
% Copyright Šimon Cieslar, Brno University of Technology, 2024-2025


%% Parse inputs
p = inputParser;
validNonNeg = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
validBinary = @(x) isnumeric(x) && isscalar(x) && any(x == [0,1]);

addRequired(p,  'chan_Signal', @(x) isnumeric(x) && ~isempty(x));
addRequired(p,  'invTestSignal', @(x) isnumeric(x) && ~isempty(x));
addRequired(p,  'fs', @(x) isnumeric(x) && isscalar(x) && (x > 0));
addParameter(p, 'SignalStart', NaN, validNonNeg);
addParameter(p, 'AutoTrim', 1, validBinary);
addParameter(p, 'ScalingFactor', 1, @(x) isnumeric(x) && isscalar(x));

parse(p, chan_Signal, invTestSignal, fs, varargin{:});
T0       = p.Results.SignalStart;
doTrim   = p.Results.AutoTrim;
scaleF   = p.Results.ScalingFactor;

%% Trim leading silence/noise from chan_Signal if requested

if doTrim
    if isnan(T0)
        % automatic detection
        threshold = 0.01 * max(chan_Signal(:).^2); % energy threshold (1% of max)
        startIdx  = detectSignalStart(chan_Signal, fs, threshold);
        fprintf('Automatic detection: Signal starts approximately at %.2f s.\n', startIdx/fs);
    else
        % manual start
        startIdx = round(T0 * fs);
        if startIdx < 1, startIdx = 1; end
        fprintf('Manual start set to %.2f s.\n', startIdx/fs);
    end
else
    startIdx = 1;
    fprintf('Trimming disabled; using full signal.\n');
end

chan_Signal = chan_Signal(startIdx:end, :);

%% Ensure that invTestSignal has the same dimensions as chan_Signal.
% if ~isequal(size(chan_Signal), size(invTestSignal))
%     chanSignalSize = size(chan_Signal); % Get the dimensions of the channel signal.
%     if isvector(invTestSignal)
%         % If invTestSignal is a vector, first ensure it is a column vector.
%         invTestSignal = reshape(invTestSignal, [], 1);
%         % Replicate the vector to match the number of elements in the remaining dimensions.
%         invTestSignal = repmat(invTestSignal, [1, prod(chanSignalSize(2:end))]);
%         % Reshape it back to the dimensions of channelSignal.
%         invTestSignal = reshape(invTestSignal, chanSignalSize);
%     else
%         % If invTestSignal is not a vector, replicate it along the missing dimensions.
%         invTestSignal = repmat(invTestSignal, [1, chanSignalSize(2:end)]);
%     end
% end
%%
% Determine the length needed for the FFT to perform linear convolution.
% For two signals of lengths nS and nInv, the linear convolution result will have length (nS + nInv - 1).
nS   = size(chan_Signal, 1);        % Number of samples in the channel signal.
nInv = size(invTestSignal, 1);      % Number of samples in the inverse test signal.
N    = nS + nInv - 1;               % Total length of the convolution result.

%%
% % Compute the FFT of both signals with the length N.
% % Specifying N as the FFT length automatically pads the signals with zeros,
% % which is necessary for linear convolution.
% fftChannel = fft(channelSignal, N);
% fftInvTest = fft(invTestSignal, N);
% 
% % Multiply the FFTs element-wise.
% % In the frequency domain, convolution corresponds to multiplication.
% fftProduct = fftChannel .* fftInvTest;
% 
% % Compute the inverse FFT to transform the result back to the time domain.
% % The result is the impulse response (IR) of the channel.
% IR = ifft(fftProduct, N);

%%
% Perform FFT-based convolution with zero-padding implicitly handled.
IR = ifft( fft(chan_Signal, N) .* fft(invTestSignal, N) );
% Scaling application
IR = scaleF * IR;

%% Helper function for signal start detection

    function startIndex = detectSignalStart(sig, fs, threshold)
    % DETECTSIGNALSTART finds the first sample index where the signal energy
    % exceeds the given threshold.
    %
    %   startIndex = DETECTSIGNALSTART(sig, fs, threshold)
    %
    %   Inputs:
    %     sig       — Input signal (vector or matrix, time along rows).
    %     fs        — Sampling frequency in Hz.
    %     threshold — Scalar energy threshold.
    %
    %   Output:
    %     startIndex — Index of the first sample above threshold (min = 1).
    
    windowSize = round(0.01 * fs);      % 100 ms window
    energy     = movmean(sig.^2, windowSize, 1);
    idx        = find(any(energy > threshold, 2), 1, 'first');
    startIndex = max(idx, 1);
    end

end