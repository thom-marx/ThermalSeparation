within ThermalSeparation.Components.Columns.BaseClasses;
record Results
  parameter Integer nSV;
  parameter Integer nSL;
  parameter Integer n;
  SI.Temperature T_l[n];
  SI.MoleFraction x_l[n,nSL];
  SI.MoleFraction x_v[n,nSV];
  SI.MoleFraction x_l_star[n,nSL];
  SI.MoleFraction x_v_star[n,nSV];
  SI.Concentration c_l[n,nSL];
  SI.Concentration c_v[n,nSV];
  SI.Pressure p_v[n + 1];
  SI.VolumeFlowRate Vdot_l[n];
  SI.VolumeFlowRate Vdot_v[n];
  Real eps_liq[n];
  Boolean startUp[n];
end Results;
