function audiodata = IR_signal_exp_sweep(dur, varargin)
% IR_SIGNAL_EXP_SWEEP generates an exponential sine‐sweep test signal and its
% inverse filter for impulse‐response measurement using the swept‐sine method.
%
%   AUDIODATA = IR_SIGNAL_EXP_SWEEP(DUR) generates an exponential (logarithmic) 
%   sine sweep signal along with its inverse filter based on the specified 
%   duration DUR (in seconds). The generated signals are intended for use 
%   in impulse response measurements following the swept-sine technique 
%   described by Farina (2000).
%
%   AUDIODATA = IR_SIGNAL_EXP_SWEEP(DUR, START_FREQ, END_FREQ) specifies the
%   sweep starting frequency (START_FREQ) and ending frequency (END_FREQ) in Hz
%   (default: 20 Hz, 20000 Hz).
%
%   AUDIODATA = IR_SIGNAL_EXP_SWEEP(DUR, START_FREQ, END_FREQ, FS) sets the
%   sampling frequency FS in Hz (default: 96000 Hz).
%
%   AUDIODATA = IR_SIGNAL_EXP_SWEEP(..., REVERSE) when REVERSE = 1 produces
%   a descending sweep (default: 0, ascending).
%
%   AUDIODATA = IR_SIGNAL_EXP_SWEEP(..., RCOS_MS) applies a raised‐cosine
%   fade‐in/fade‐out of RCOS_MS (in milliseconds) (default: 15 ms) to smooth 
%   the sweep edges.
%
%   AUDIODATA = IR_SIGNAL_EXP_SWEEP(..., DOPLOT) if DOPLOT = 1, displays the
%   time‐domain waveforms and spectrograms via plotIRVisualization (default: 0).
%
%   Reference:
%       Farina, A. (2000). Simultaneous Measurement of Impulse Response and 
%       Distortion With a Swept-Sine Technique. AES Preprint.
%       [Available at: https://www.researchgate.net/publication/2456363_Simultahneous_Measurement_of_Impulse_Response_and_Distortion_With_a_Swept-Sine_Technique]
%
% Copyright Šimon Cieslar, Brno University of Technology, 2024-2025

%% Parse Input Arguments
[dur, start_freq, end_freq, fs, reverse, rcos_ms, doPlot] = parseInputs(dur, varargin{:});

%% Setup Basic Parameters
SampleInterval = 1/fs;        % Sampling interval in seconds.
amplitude      = 0.5;         % Base amplitude of the sweep signal.
scale_inv      = 1;           % Enable allpass compensation.

%% Generate Exponential Sweep Signal
% According to Farina (2000, p.6), the exponential sweep is defined by two constants K and L.
%
% The sweep phase is given by:
%     phi(t) = K * [exp(L * t) - 1]
%
% For a sweep starting at ω1 (2π·start_freq) and ending at ω2 (2π·end_freq) over a duration T,
% the constants K and L are derived (see Farina, 2000, p.6) as:
%
%     K = (T * ω1) / ln(ω2/ω1)
%     L = ln(ω2/ω1) / T
%
% In the code, T is "dur" and the angular frequencies are computed as:
w1 = 2 * pi * start_freq;  % Starting angular frequency (rad/s)
w2 = 2 * pi * end_freq;    % Ending   angular frequency (rad/s)

% Compute K and L (cf. Farina, 2000, p.6)
K = (dur * w1) / log(w2 / w1);  
L = log(w2 / w1) / dur;         

% Create time vector for the signal duration.
nSamples = dur / SampleInterval;
t        = linspace(0, dur - SampleInterval, nSamples);

% Compute the phase function of the sweep signal.
phi = K * (exp(t * L) - 1);  % Exponential phase function (Farina, 2000, p.6)

% Compute instantaneous frequency (can be used for amplitude shaping).
freq = K * L * exp(t * L);

% Calculate an amplitude envelope to compensate for the energy variation over frequency.
% This amplitude shaping (≈6 dB/octave reduction) is introduced to equalize the energy,
% as described in Farina (2000, p.6-7).
amp_env = 10.^((log10(amplitude)) * log2(freq / freq(1)));

% Generate the exponential sweep signal using the sine function.
S = amplitude * sin(phi);

%% Apply Raised Cosine Window (Fade-In and Fade-Out)
% To smooth the edges of the signal and avoid abrupt transitions, apply a raised cosine (Hann) window.
% Compute the number of samples corresponding to the fade duration (rcos_ms) based on the sampling frequency.
rcos_len = (rcos_ms * 1e-3) * fs;
sig_len = length(S);

% Ensure that the combined fade durations (fade-in and fade-out) do not exceed the signal length.
if 2 * rcos_len > sig_len
    warning(['Specified raised cosine duration is too high relative ' ...
        'to signal length. Reducing fade duration to 1 % of the signal length.']);
    rcos_len = floor(0.01 * sig_len);
end

% Generate the raised cosine window (using a Hann window).
rcoswin = hann(2 * rcos_len).';

% Apply the window to the beginning and end of the signal.
S = [S(1:rcos_len) .* rcoswin(1:rcos_len), ...
     S(rcos_len+1:sig_len-rcos_len), ...
     S(sig_len-rcos_len+1:sig_len) .* rcoswin(rcos_len+1:end)];

%% Generate Inverse Filter
% The inverse filter is obtained by time-reversing S and applying the amplitude envelope.
% This procedure (time reversal and amplitude correction) is detailed in Farina (2000, p.6-7).
Sinv = fliplr(S) .* amp_env;

%% Correct for Allpass Delay in the Inverse Filter
% To compensate for the phase delay introduced by the sweep (which acts as an allpass filter),
% the inverse filter is corrected via an FFT-based phase compensation.
t_sinvfft = linspace(0, sig_len-1, sig_len);
Sinvfft   = fft(Sinv);
Sinvfft   = Sinvfft .* exp(1i * 2*pi * t_sinvfft * (sig_len-1)/sig_len);
Sinv      = real(ifft(Sinvfft));

%% Reverse Signal if Requested
% If the REVERSE flag is set, reverse both the sweep and inverse filter.

if scale_inv == 1
   fftS      = fft(S);
   mid_freq  = (start_freq + end_freq)/2;
   index     = round(mid_freq/(fs/sig_len));
   const1    = abs(conj(fftS(index))/(abs(fftS(index))^2));
   const2    = abs(Sinvfft(index));
   % Scaling factor is applied in convolution
   IRscalingfactor = const1/const2;
else
    IRscalingfactor = 1;
end

if reverse
    S    = fliplr(S);
    Sinv = fliplr(Sinv);
    % Reversal is discussed in Farina (2000, p.6)
end

%% Normalize the RMS of the final singal
targetRMS = 0.07;  % Target RMS level for the active signal
S    = S    * (targetRMS / rms(S));
Sinv = Sinv * (targetRMS / rms(Sinv));

%% Assemble Output Structure
audiodata.audio  = S.';
audiodata.audio2 = Sinv.';
audiodata.fs     = fs;
audiodata.inarg  = {dur, start_freq, end_freq, fs, reverse, rcos_ms, doPlot};
audiodata.IRscalingfactor = IRscalingfactor;

%% Optional Visualizationha
    if doPlot
        titleStr = sprintf('ESS %ds %d–%dHz at %dHz', dur, start_freq, end_freq, fs);
        plotIRVisualization(audiodata, titleStr);
    end
end

%% Helper Function: parseInputs
function [dur, start_freq, end_freq, fs, reverse, rcos_ms, doPlot] = parseInputs(dur, varargin)
% PARSEINPUTS Parses input arguments for IR_signal_exp_sweep.
%
%   [DUR, START_FREQ, END_FREQ, FS, REVERSE, RCOS_MS, DOPLOT] = parseInputs(DUR, ...)
%
%   Required:
%       DUR         - Duration of the signal in seconds.
%
%   Optional (in order):
%       START_FREQ  - Starting frequency in Hz (default: 20 Hz).
%       END_FREQ    - Ending frequency in Hz (default: 20000 Hz).
%       FS          - Sampling frequency in Hz (default: 96000 Hz).
%       REVERSE     - 0 for ascending sweep, 1 for descending sweep (default: 0).
%       RCOS_MS     - Duration for fade-in and fade-out in milliseconds (default: 15 ms).
%       DOPLOT      - 0 = no plots, 1 = visualize via plotIRVisualization (default: 0).
%
    p = inputParser;
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
    validFlag         = @(x) isnumeric(x)&&isscalar(x)&&ismember(x,[0,1]);

    addRequired(p, 'dur', validScalarPosNum);
    addOptional(p, 'start_freq', 20, validScalarPosNum);
    addOptional(p, 'end_freq', 20000, validScalarPosNum);
    addOptional(p, 'fs', 96000, validScalarPosNum);
    addOptional(p, 'reverse', 0, validFlag);
    addOptional(p, 'rcos_ms', 15, validScalarPosNum);
    addOptional(p, 'doPlot', 0, validFlag);

    parse(p, dur, varargin{:});
    
    dur        = p.Results.dur;
    start_freq = p.Results.start_freq;
    end_freq   = p.Results.end_freq;
    fs         = p.Results.fs;
    reverse    = p.Results.reverse;
    rcos_ms    = p.Results.rcos_ms;
    doPlot     = p.Results.doPlot;
end