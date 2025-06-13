# Extremum Seeking Control for Resonant Frequency Seeking
Extremum seeking control (ESC) is implemented to obtain the resonant frequency of an underdamped mass-spring-damper system. For this purpose, a control architecture is designed for ESC to modulate the frequency of the input force and measure an estimate of the amplitude of the mass displacement. 
Continuous-time and discrete-time implementations are shown
A vanishing dither algorithm for converging ESC output (control input) is featured in the discrete-time implementation.
Documentation of the algorithms implemented can be found in the *Extremum_Seeking_for_Resonant_Frequency_Seeking.pdf* file.

## Main Files
* **msd_resonant_frequency_seeking_CT_run.m** : Program for setting up and running the continuous-time, extremum seeking control for resonant frequency seeking simulation. Modulation dither amplitude is arbitrarily made zero after some time.
* **msd_resonant_frequency_seeking_CT.slx** : Simulink file for closed-loop simulation of the mass-spring-damper system and the continuous-time, extremum seeking control for resonant frequency seeking algorithm.
* **msd_resonant_frequency_seeking_DT_run.m** : Program for setting up and running the discrete-time, extremum seeking control for resonant frequency seeking simulation. Modulation dither amplitude is arbitrarily made zero after some time.
* **msd_resonant_frequency_seeking_DT.slx** : Simulink file for closed-loop simulation of the mass-spring-damper system and the discrete-time, extremum seeking control for resonant frequency seeking algorithm.
* **msd_resonant_frequency_seeking_DTvanishing_dither_run.m** : Program for setting up and running the discrete-time, extremum seeking control for resonant frequency seeking simulation, with a vanishing modulation dither.
* **msd_resonant_frequency_seeking_DTvanishing_dither.slx** : Simulink filefor closed-loop simulation of the mass-spring-damper system and the discrete-time, extremum seeking control for resonant frequency seeking algorithm, with a vanishing modulation dither.

## Support Files
* **Extremum_Seeking_for_Resonant_Frequency_Seeking.pdf** : Documentation of the continuous-time and discrete-time Extremum Seeking Control for Resonant Frequency Seeking algorithms.
* **msd_frequency_response_magnitude.m** : Plots the magnitude of the frequency response of the mass-spring-damper system.
* **nonlinear_oscillator_frequency_switch_test.m** : Tests the continous-time and discrete-time nonlinear oscillators used for frequency modulation in the Resonant Frequency Seeking simulations.
