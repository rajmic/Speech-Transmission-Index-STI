# Speech-Transmission-Index-STI-
Matlab implementation of STI (Full STI, STIPA, direct and indirect approaches, correction factors)

**Speech Transmission Index (STI)** ([STI](https://en.wikipedia.org/wiki/Speech_transmission_index))
is a metric ranging between 0 and 1 predicting the speech intelligibility when speech is passed through a transmission channel, defined in the [IEC&nbsp;60268-16](https://webstore.iec.ch/publication/26771) standard [[1]](#1). It is based on an analysis of the amplitude modulations, which simulate speech signals.

The quality of speech transmission and the likelihood of intelligibility of syllables, words, and sentences being comprehended for native speakers with healthy hearing can be represented by the following table.

| STI value | Quality according to<br>IEC 60268-16 | Intelligibility<br>of syllables in % | Intelligibility<br>of words in % | Intelligibility<br>of sentences in % |
|:----------------:|:---------:|:-------------:|:-------------:|:--------------:|
| 0 &ndash; 0.3    | bad       | 0 &ndash; 34  | 0 &ndash; 67  | 0 &ndash; 89   |
| 0.3 &ndash; 0.45 | poor      | 34 &ndash; 48 | 67 &ndash; 78 | 89 &ndash; 92  |
| 0.45 &ndash; 0.6 | fair      | 48 &ndash; 67 | 78 &ndash; 87 | 92 &ndash; 95  |
| 0.6 &ndash; 0.75 | good      | 67 &ndash; 90 | 87 &ndash; 94 | 95 &ndash; 96  |
| 0.75 &ndash; 1   | excellent | 90 &ndash; 96 | 94 &ndash; 96 | 96 &ndash; 100 |

This MATLAB project provides an open-source implementation of STI methods based on an implementation developed in the master’s thesis "_[Computing Speech Transmission Index in Matlab](https://github.com/Cieslar-Simon/STI)_" (Š. Cieslar, Brno University of Technology, 2025) and on our earlier Matlab implementation of STIPA [[2]](#2)

**Full STI vs. STIPA**: The full STI procedure uses 98 separate test signals covering 7 octave frequency bands and 14 modulation frequencies, requiring roughly a 15-minute measurement. **STIPA** (Speech Transmission Index for Public Address systems) is a simplified, standardized variant of STI that uses a single composite test signal with only two modulation frequencies per band, reducing measurement time to about 15–25 seconds. STIPA provides a quicker estimation of STI with minimal accuracy loss, making it practical for field measurements. This repository implements both the **direct STI methods** – including the exhaustive Full STI and the faster STIPA – and the **indirect method** based on impulse response analysis. The _direct_ method involves passing a known modulated noise signal through the system and analyzing the recorded output, while the _indirect_ method computes STI by analyzing the system’s measured impulse response (IR). (The indirect IR-based method assumes a linear time-invariant system; it may be less accurate if the audio path includes non-linear elements like companders or limiters.) By providing implementations of all these approaches, the project allows users to compare results and understand how different measurement techniques can be used to assess speech intelligibility.

## Features
+ **Full Direct STI Calculation**: Generation of the full STI test signal (7 octave bands × 14 modulations) and analysis of recorded signals to compute the STI. This includes calculating the Modulation Transfer Function (MTF) for each band/modulation frequency and applying standard STI weighting and redundancy corrections. The code can produce the long-format test signal and analyze its playback recording to yield an STI score along with intermediate metrics (modulation indices, per-band intelligibility scores, etc.).
+ **STIPA Method**: Implementation of the **STIPA** simplified test. The repository can generate the standard STIPA test signal (modulated noise with two modulation frequencies per band) and compute STI from a recorded STIPA test playback. This provides a fast estimate of intelligibility using the widely adopted STIPA procedure (per IEC 60268-16). It includes functions for creating the STIPA signal of a desired length and sampling rate, as well as analyzing the recorded response to output the STI value and modulation reduction per octave band.
+ **Indirect STI from Impulse Response**: Tools to calculate STI based on an IR analysis. Given a measured IR of the environment or audio chain, the code computes the modulation transfer function by processing the IR (for example, by deriving the frequency response and applying octave-band modulation analysis). This yields an STI prediction without needing the direct playback of modulated noise. An exponential sine sweep signal generation (<code>IR_signal_exp_sweep.m</code>) is provided to facilitate IR measurements, and the measured IR can be fed into the <code>sti_ir.m</code> function to compute STI.
+ **Signal Generation and Processing**: Functions to generate all necessary test signals and to process recorded signals:
  + <code>generateFullSTISignal.m</code> creates the full STI multi-band test signal according to the standard (all modulation frequencies in all octave bands).
  + <code>generateStipaSignal.m</code> produces a STIPA test signal of specified duration and sample rate, containing the proper modulation patterns [github.com](https://github.com/zawi01/stipa#:~:text=The%20full%20STI%20model%20consists,15%20s%20to%2025%20s).
  + <code>IR_signal_exp_sweep.m</code> generates an exponential swept-sine signal for measuring impulse responses (a common method to capture IR while minimizing distortion).
  + Utility functions like <code>convolution.m</code> (for convolving signals with impulse responses or filters) are included to assist in signal processing tasks.
+ **Octave Band Analysis**: The implementation follows the standard’s use of 7 octave frequency bands (centered at 125&nbsp;Hz up to 8&nbsp;kHz). The function <code>octave_band_analysis.m</code> filters signals or IRs into these octave bands and computes band-specific metrics. This is used both in direct and indirect computations to evaluate modulation depth reduction per band in the presence of noise and reverberation.
+ **Comprehensive Analysis Scripts**: Three demonstration scripts are provided to show end-to-end usage:
  + <code>demonstration_fullsti.m</code> - guides the user through generating the full STI signal, playing or simulating its transmission, analyzing the received signal with <code>fullsti.m</code>, and obtaining the STI result. It also visualizes the modulation transfer matrix or intermediate results.
  + <code>demonstration_stipa.m</code> - demonstrates the STIPA procedure: it can generate or load a STIPA test signal, help simulate or record its playback, then use <code>stipa.m</code> to compute the STI. This script may plot the octave-band modulation indices and display the final STI value in a table or console.
  + <code>demonstration_sti_ir.m</code> - shows how to use a measured impulse response to compute STI via the indirect method. For example, it might load an example IR (provided from a real room measurement) or use a synthetic IR, and then call <code>sti_ir.m</code> to calculate the STI. The script can then compare this indirect result with direct measurements (if available) for validation.
+ **Visualization and Output**: The project includes functions to help interpret the results:
  + <code>displayPlotSTI.m</code> generates plots of the modulation reduction matrix or bar charts of the transmission index per band, helping users visualize how each octave band contributes to the overall STI.
  + <code>displayTableSTI.m</code>prints or returns a formatted table of STI results, possibly including per-octave-band STI values and overall STI (and maybe the CIS (common intelligibility scale) or quality rating).
  + <code>plotSignalWaveform.m</code> and <code>plotIRVisualization.m</code> are utility functions to plot time-domain waveforms of test signals or impulse responses and possibly their envelopes or frequency responses. These help users verify the signals and IRs (for example, checking the quality of a recorded test signal or the decay of an impulse response).
+ **Calibration Support**: To test the demonstration scripts, a reference calibration tone file <code>calibrator94dB.wav</code> is included, containing a 1 kHz tone at 94&nbsp;dB SPL. This can be used to calibrate the recording system so that the recorded signal levels correspond to real sound pressure levels. By recording this tone with the measurement microphone and adjusting gain such that the analysis recognizes it as 94 dB, subsequent STI measurements can be referenced to absolute SPL. Proper calibration ensures that noise levels and speech levels are correctly accounted for in the STI calculation (per the standard’s requirements).

## Structure of the Repository
The repository is organized into MATLAB scripts, functions, and data files, with a dedicated folder for measurement data:
+ **Root Directory**: Contains the main analysis functions, scripts, and supporting utilities:
  + **STI computation functions**: <code>fullsti.m</code> (computes Full STI from a recorded signal), <code>stipa.m</code> (computes STI from a STIPA test signal recording), and <code>sti_ir.m</code> (computes STI from an impulse response). These are the core algorithms implementing the STI calculations for each method.
  + **Signal generation**: <code>generateFullSTISignal.m</code> and <code>generateStipaSignal.m</code> create the test signals for Full STI and STIPA respectively. These functions follow the standard modulation specifications to output a waveform that can be played over a loudspeaker for measurement or used in simulations.
  + **Demonstration scripts**: <code>demonstration_fullsti.m</code>, <code>demonstration_stipa.m</code>, <code>demonstration_sti_ir.m</code> illustrate how to use the above functions. Users can run these scripts in MATLAB to see example workflows (generating signals, reading measurement files, computing STI, and displaying results). Each script is documented with steps for clarity.
  + **Utilities**:
    + <code>octave_band_analysis.m</code> - filters an input signal or IR into the 7 octave bands and computes relevant band-level metrics (e.g., RMS levels or modulation indices in each band).
    + <code>addSilenceGaps.m</code> - inserts silent intervals into a signal (e.g., between bursts of noise). This can be useful for formatting test signals or separating sections in a composite signal.
    + <code>convolution.m</code> - performs convolution of two signals (used for applying an impulse response to a test signal, if simulating a transmission through an acoustic channel within MATLAB).
    + Plotting and display helpers: <code>displayPlotSTI.m</code>, <code>displayTableSTI.m</code> (for results visualization), <code>plotIRVisualization.m</code> (for plotting impulse response characteristics), and <code>plotSignalWaveform.m</code> (for plotting audio waveforms).
  + **Test data files**: Example data is provided to facilitate quick testing:
    +  <code>stipa_test_signal.wav</code> - an example STIPA test signal (generated noise with proper modulation). This can be used for quick experiments or played back for measurement. Users can also generate a fresh one using <code>generateStipaSignal.m</code>.
    +  <code>stipaMeasurement.wav</code> - an example recorded STIPA signal captured during the measurement campaign (playing the STIPA test signal through a system and recording at a listener position). This file allows users to try out the <code>stipa.m</code> analysis without needing to immediately perform their own recording.
    +  <code>calibrator94dB.wav</code> - a calibration tone (1 kHz sine at 94 dB SPL). This short recording is used to calibrate measurement equipment or verify that the analysis correctly interprets levels. Playing and recording this file in a measurement setup should register as 94 dB in the analysis, ensuring correct level calibration for STI calculations.

## Usage Instructions

**Prerequisites**: To use this repository, you need MATLAB (the code was developed and tested in MATLAB R2024b). No specialized toolboxes are strictly required for basic functionality – the code uses standard MATLAB functions for signal processing (octave-band filtering is implemented manually in <code>ctave_band_analysis.m</code>). However, having the Signal Processing Toolbox or Audio Toolbox may be beneficial for audio I/O and visualization but is not mandatory. Ensure your MATLAB path includes the repository files (you can achieve this by running <code>addpath(genpath('<path-to-repo>'))</code> or by opening MATLAB in the repository folder).

**Running a Basic STIPA Analysis**: A simple way to get started is to run the provided STIPA demonstration:
1. Launch MATLAB and navigate to the repository folder.
2. Open and run <code>demonstration_stipa.m</code>. This script will:
  + Generate a STIPA test signal (or load the provided <code>stipa_test_signal.wav</code>).
  + If using provided data: it may load the example recorded response <code>stipaMeasurement.wav</code> to simulate a real measurement.
  + Call the <code>stipa()</code> function to compute the Speech Transmission Index from the recorded signal.
  + Display the resulting STI value, and plot the modulation depth loss per octave band using <code>displayPlotSTI.m</code> or output a table via <code>displayTableSTI.m</code>.
  + The script will output a final STI score in the MATLAB console and possibly open figures illustrating the intelligibility per frequency band.
3. Review the console output and plots. For example, you might see a plot of 14 modulation indices across 7 bands and a printed table of STI contributions per band, culminating in the overall STI rating (e.g., _STI = 0.65 (Good)_).

**Performing a Full STI Measurement**: To use the full STI method:
1. Use <code>generateFullSTISignal.m</code> to create the full STI test signal. By default, it may generate the standard 15-minute signal. You can specify parameters such as sampling rate or duration if needed.
2. Play this signal through the system or environment under test (e.g., through a loudspeaker in a room). Make sure to record the output with a suitable microphone or recording device.
3. Save the recorded audio (it should contain all 98 modulation combinations in sequence). Due to length, ensure no interruptions occur and the recording device can handle the duration.
4. Use <code>fullsti.m</code> to analyze the recorded waveform. For example: <code>[STI_value, details] = fullsti('recordedFullSTI.wav', fs);</code> (where fs is the sampling rate). This will apply octave filtering, compute modulation transfer ratios for each frequency/modulation combination, and calculate the STI.
5. Review the output. The function may return the overall STI and optionally a structure or matrix of modulation indices. You can use <code>displayPlotSTI.m</code> to visualize the full STI modulation matrix or <code>displayTableSTI.m</code> to see a summary of results.

**Using the Indirect IR Method**: If you have an impulse response of the system:
1. Ensure the impulse response is available as a WAV or MATLAB array.
2. Load the impulse response in MATLAB (e.g., <code>[ir, fs] = audioread('RIR.wav');</code>).
3. Run the <code>sti_ir.m</code> function with this IR: <code>STI_val = sti_ir(ir, fs);</code>. The function will likely internally filter the IR into octave bands, compute the modulation frequency response (either by convolving the IR with modulated noise or using the IR’s frequency response to derive MTF according to the standard’s prescribed method), and then compute the STI.
4. The result <code>STI_val</code> is the estimated speech transmission index. You can compare this with any direct measurement at the same position. The demonstration script <code>demonstration_sti_ir.m</code> can be used as a guide; it may automatically load one of the provided IRs and perform the above steps, then possibly compare the indirect STI to a STIPA result at that position.
5. If you need to measure an IR yourself, you can use <code>IR_signal_exp_sweep.m</code> to generate a sweep signal, play it in the environment, record the response, then deconvolve (the script or repository might have instructions for deconvolving the recorded sweep to get the IR).

**Calibration (Optional)**: If absolute accuracy is required (especially in presence of noise), use the calibration tone:
1. Play the <code>calibrator94dB.wav</code> file through your playback system with a standard 94 dB SPL calibrator (or use it to adjust your system output).
2. Record it with the microphone in the measurement chain.
3. Adjust your recording gain such that analyzing this recording yields ~94 dB SPL in the internal calculations. This might involve ensuring that the waveform amplitude corresponds to the calibration level expected by the STI algorithm (some STI calculations need the absolute speech and noise levels).
4. Once calibrated, proceed with STI measurements. This step ensures that the “speech level” and “noise level” inputs to STI calculation are accurate in dB SPL, as required by the IEC standard. (If calibration is not performed, the STI computation will still work, but the results are relative – which is usually acceptable as long as the test signal playback level is within the standard’s specified range.)

**General Notes**: All analysis functions output or plot results that align with the STI standard’s outputs. For example, you can expect the STI score and possibly intermediate measures like _Modulation Index (m)_ for each octave and modulation frequency, _MTF_(modulation transfer function values), and signal-to-noise ratios if noise is involved. The code is documented with comments to help understand the steps. If you wish to modify the code (e.g., to use different noise signals or to test only certain octave bands), the comments and modular structure will assist in making targeted changes.

## References
<a id="1">[1]</a> 
International Electrotechnical Commission,“Sound system equipment – Part 16: Objective rating of speech intelligibility by speech transmission index,” 2020. Number IEC 60268-16:2020, edition 5.0.

<a id="2">[2]</a> 
Záviška, P., Rajmic, P., Schimmel, J. _MATLAB Implementation of STIPA_. AES Europe Conference, 2024 – An open-source project implementing the STIPA method for speech intelligibility testing.

&copy; Šimon Cieslar, Pavel Záviška, Jiří Schimmel, Pavel Rajmic, Brno University of Technology, 2023&ndash;2025
