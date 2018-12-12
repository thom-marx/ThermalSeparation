within ThermalSeparation.Components.GasLiquidVolumes;
model LSF_noDeltaP "lean solvent flash where outlet pressure = inlet pressure"
extends ThermalSeparation.Components.GasLiquidVolumes.BaseClasses.BaseLSF;
equation
p=p_out;
initial equation

end LSF_noDeltaP;
