within ThermalSeparation.Components.LiquidVolumes;
record ResultsLiquidVolumes
parameter Integer nSL;
SI.MoleFraction x_l[nSL];
SI.Height level;
SI.Volume V_liq;
SI.Concentration c_l[nSL];
SI.Temperature T;
end ResultsLiquidVolumes;
