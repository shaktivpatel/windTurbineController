# windTurbineController

Three different pitch angle controllers (PID observer & state-feedback) for a linearized wind turbine model

Objectives
---------

* To model a linearized wind turbine using MATLAB/Simulink
* To design three different pitch angle controllers: PID, state-feedback and observer-based.
* To conduct stability analysis, simulate dynamic system responses, and visualize controller performance.

Background
--------

<img align="right" width="266" height="145" alt="image" src="https://github.com/user-attachments/assets/2e8c104e-b5bf-4ef2-9712-29522c8175de" />
Based on the different operating regions of a wind turbine, in region 3 a controller is needed to ensure that the generated power is capped at the rated power. The pitch angle controller that is designed to monitor and adjust the wind turbine's rotor blade angle which directly affects the power generation.


Results
------
The follow inputs were tested with the designed controllers:
* Randomized wind speed profile (4.7m/s-8.1m/s)
* Sinusoidal input
* Step input
* Impulse input

