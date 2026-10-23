# Virtual Wind Chime

An audio signal processing study exploring how irregular excitation and decaying resonances can create an evolving wind-chime soundscape. The repository combines the original MATLAB analysis with a playable browser interpretation for the portfolio.

## Project credits

Final project for the Audio Signal Processing course in the master's programme in AAT at Aalto University. Original team: Adrain Kalazi, Maria Pikhtina and Sara Ghorbani. Advisor: Dr. Vesa Välimäki. Names follow the original repository credits.

## Play the prototype

Open `docs/index.html` in a current browser and press **Start listening**. Alternatively, serve the repository locally:

```sh
python3 -m http.server 8000
```

Then open `http://localhost:8000/docs/`. No package installation or build step is required. Fonts are fetched from Google Fonts; system fonts provide a fallback.

Click or tap any tube, or press keys **1–5**, to strike it. **Wind intensity** changes the stochastic strike rate and strength. **Send a gust** creates a short burst of activity. **Resonance** scales modal decay times without changing pitch; **Output** controls listening level. The mode inspector displays the five source frequencies for the selected tube. Headphones reveal the stereo placement. Pause fades active voices; switching away from the page stops playback.

## Sound design

The browser renderer sums five damped sinusoidal modes per strike. It uses the hard-coded `tube_modes`, `gains` and `pole_radii` from `windchimes.m`. These frequency sets are not forced into a harmonic series or a musical scale. The unusual 64.6 Hz entry for Tube 1 is retained as source data, not validated as a measured resonance.

| Tube | Mode 1 (Hz) | Mode 2 (Hz) | Mode 3 (Hz) | Mode 4 (Hz) | Mode 5 (Hz) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 64.6 | 1055.1 | 2659.4 | 3617.6 | 4974.2 |
| 2 | 592.16 | 1571.9 | 2960.8 | 4199 | 5437.1 |
| 3 | 592.16 | 1571.9 | 3456.1 | 5351 | 6632.2 |
| 4 | 592.16 | 1571.9 | 2896.2 | 3886.7 | 5254.1 |
| 5 | 796.73 | 1755 | 2734.7 | 4252.8 | 5221.8 |

Relative modal gains are `[0.0550, 0.0877, 0.0719, 0.0586, 0.0331]`. Pole radii are `[0.9998, 0.9997, 0.9992, 0.9995, 0.9980]`, applied across all tubes. Mode height in the inspector represents these gains, not a measured or live spectrum.

For each mode, the browser envelope uses a 2 ms onset followed by exponential decay. The decay constant is `tau = -1 / (44100 * ln(r))`, scaled by the resonance control. A nominal 44.1 kHz reference rate is an explicit assumption: the source recording's actual sample rate is unavailable. The browser runs at the device AudioContext sample rate; modes at or above 45% of that rate are omitted. Voices fade to zero after approximately 60 dB of decay.

Each tube has a fixed stereo position. A master compressor controls accumulation from overlapping strikes. No reverberation or wind-noise layer is added, keeping modal decay audible. Wind is represented by irregular event timing and variable excitation strength. The scheduler, stereo mapping, 2 ms onset and motion graphics are prototype additions, not recovered original functions. Geometry and motion are illustrative rather than a calibrated mechanical model.

## Original MATLAB workflow and limitations

The script reads the opening two seconds of `wood.wav`, computes a 4096-point spectrogram with 90% overlap, and attempts modal-frequency and gain estimation across five manually selected time regions. It then specifies a 25-second simulation and calls wind, energy and clapper helper functions.

The repository does not contain `wood.wav`, `simulate_wind_force`, `update_system_energy`, `simulate_clapper_movement` or `synthesize_tube_sound`. The final synthesis call is commented out. The MATLAB script therefore cannot reproduce the complete audio pipeline from this repository alone, and the browser prototype is not a bit-exact port or a validation of the original results.

Several analysis details require correction before using the MATLAB estimator for quantitative work. `10*log10(abs(S))` is neither standard amplitude dB nor power dB; use `20*log10(abs(S))` for magnitude or `10*log10(P)` for power with an explicit reference. Summing logarithmic spectra is not linear power integration. `findpeaks` distance is in samples of the supplied spectral vector, so a desired separation in Hz must be converted to frequency bins. The script takes the first detected peaks, not necessarily the strongest modes.

The parabola expression `-p(2)/(2*p(1))` returns the peak's frequency offset, not its amplitude. The browser uses the separate hard-coded gains instead. For a pole-radius amplitude envelope, `r = 10^(-3/(T60*Fs))` gives 60 dB amplitude decay in T60 seconds. The script's example exponent of -6 gives 120 dB amplitude decay over that duration. The hard-coded pole radii and the example T60 array are separate parameter sets.

## Related research

Teemu Lukkari and Vesa Välimäki (2004). **Modal Synthesis of Wind Chime Sounds with Stochastic Event Triggering.** Proceedings of the 6th Nordic Signal Processing Symposium (NORSIG 2004), pp. 208–211.

[Aalto University publication record](https://research.aalto.fi/en/publications/modal-synthesis-of-wind-chime-sounds-with-stochastic-event-trigge/)

This is a closely relevant research reference. The original README and script do not identify the paper used in the course project, so its status as the project's original reference remains unconfirmed. The browser implementation should not be described as an exact implementation of this paper.

## Verification

JavaScript syntax was checked with `node --check docs/app.js`. Browser layout, interaction and audio playback have not yet been verified in the development environment because its preview service was unavailable. Test desktop and mobile playback before using the demo as a finished portfolio release.
