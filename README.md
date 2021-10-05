# AverageCSI_EnergyBeamforming
This compilation of MatLab scripts is used to generate the results illustrated in [REF].
The scripts description is as follows.

full-CSI WET precoding:
  - CSIBeamf_SDP (solution of P2 eq.4 in [REF])

average-CSI WET precoding:
  - CSIBeamf_SDP (solution of P2 eq.4 in [REF]) but using the average CSI information only
  - avgCSI_MRT (MRT precoder obtained from eq.6 and the power allocation from Algorithm 1 in [REF])


Scripts for figures generation:

  - ToySimulation_LOS_Factor.m (It allows generating a figure similar to [REF, Fig.2]: Average worst-case RF energy available at the user devices) 

  - ToySimulation_Number_Antennas.m (It allows generating a figure similar to [REF, Fig.3]: Average worst-case RF energy available at the user devices as a function of the number of PB's antennas)

  - ToySimulation_Rotation_Angle.m (It allows generating a figure similar to [REF, Fig.4]: Average worst-case RF energy available at the user devices as a function of the PB's rotation angle) 


References:

[REF] - O. L. A. López, F. A. Monteiro, H. Alves, R. Zhang and M. Latva-aho, "A Low-Complexity Beamforming Design for Multiuser Wireless Energy Transfer," in IEEE Wireless Communications Letters, doi: 10.1109/LWC.2020.3020576.
