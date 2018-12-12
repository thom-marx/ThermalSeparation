within ThermalSeparation.PhaseEquilibrium;
partial model BasePhaseEquilibrium
  parameter Integer nS "number of substances" annotation(Dialog(enable=false));
  constant Integer nSL= MediumLiquid.nSubstance "number of substances";
  constant Integer nSV= MediumVapour.nSubstance "number of substances";
  parameter Real factor_K[nSL] = fill(1,nSL)
    "calibration factor for phase-equilibrium";
  parameter Integer mapping[nS,2]
    "parameter to map the different medium vectors one to another"  annotation(Dialog(enable=false));
  replaceable package MediumVapour =
      ThermalSeparation.Media.BaseMediumVapour  annotation(Dialog(enable=false));
  replaceable package MediumLiquid =
      ThermalSeparation.Media.BaseMediumLiquid  annotation(Dialog(enable=false));
  input SI.Temperature T;
  input SI.MoleFraction x_l[nSL];
  input SI.MoleFraction x_v[nSV];
  input SI.MoleFraction x_vap_liq[nS];
  input SI.MolarVolume v_v;
  input SI.Pressure p;
  input SI.Pressure p_sat[nSL];
  output Real K[nS];
  output SI.Pressure p_bubble "bubble pressure - used for startup";

end BasePhaseEquilibrium;
