within ThermalSeparation.PressureLoss.PlateColumn;
model ConstantVolumeFlowRate
  "volume flow rate is constant and equal to the flow rate at the inlet - not useful for rectification"
  extends BasicPressureLossPlate;
   // output SI.Pressure dp_1;
  // input SI.VolumeFlowRate Vdot_in;
equation

for j in 1:n loop
  if empty[j] then
  Vdot[j] = 0;
  else
  Vdot[j] = Vdot_in;
  end if;
end for;

p_v_in = p[1];

end ConstantVolumeFlowRate;
