function signal = generateSweptSineSignal(duration, varargin)
% GENERATESWEPTSINESIGNAL generates an exponential sine‐sweep test signal
% using IR_SIGNAL_EXP_SWEEP function for impulse‐response measurement
%
%   SWEPTSINE_SIGNAL = GENERATESWEPTSINESIGNAL(DURATION) generates an 
%   exponential (logarithmic) sine sweep signal based on the specified 
%   duration DURATION (in seconds). The generated signals are intended for  
%   use in impulse response measurements following the swept-sine technique 
%   described by Farina (2000).
%
%   SWEPTSINE_SIGNAL = GENERATESWEPTSINESIGNAL(DURATION, START_FREQ, END_FREQ) 
%   specifies the sweep starting frequency (START_FREQ) and ending frequency 
%   (END_FREQ) in Hz (default: 20 Hz, 20000 Hz).
%
%   SWEPTSINE_SIGNAL = GENERATESWEPTSINESIGNAL(DURATION, START_FREQ, END_FREQ, FS) 
%   sets the sampling frequency FS in Hz (default: 96000 Hz).
%
%   SWEPTSINE_SIGNAL = GENERATESWEPTSINESIGNAL(..., REVERSE) when 
%   REVERSE = 1 produces a descending sweep (default: 0, ascending).
%
%   SWEPTSINE_SIGNAL = GENERATESWEPTSINESIGNAL(..., RCOS_MS) applies a 
%   raised‐cosine fade‐in/fade‐out of RCOS_MS (in milliseconds) (default: 
%   15 ms) to smooth the sweep edges.
%
%   SWEPTSINE_SIGNAL = GENERATESWEPTSINESIGNAL(..., DOPLOT) if DOPLOT = 1, 
%   displays the time‐domain waveforms and spectrograms via 
%   plotIRVisualization (default: 0).
%

    audiodata = IR_signal_exp_sweep(duration, varargin{:});
    signal = audiodata.audio;
end