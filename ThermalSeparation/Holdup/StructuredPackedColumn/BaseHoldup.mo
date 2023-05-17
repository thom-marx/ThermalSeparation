within ThermalSeparation.Holdup.StructuredPackedColumn;
partial model BaseHoldup
  parameter Integer n(min=1) annotation(Dialog(enable=false));

    replaceable record Geometry =
      ThermalSeparation.Geometry.StructuredPackedColumn.Geometry                constrainedby ThermalSeparation.Geometry.StructuredPackedColumn.Geometry;

  replaceable package MediumLiquid = Media.BaseMediumLiquid;
  replaceable package MediumVapour = Media.BaseMediumVapour;
  input MediumLiquid.ThermodynamicProperties propsLiq[n];
  input MediumVapour.ThermodynamicProperties propsVap[n];

  input SI.Pressure p_v[n+1];

  input Real eps_liq[n];
  output Real hu_stat[n];
  output SI.VolumeFlowRate Vdot[n](start=fill(0.03,n));
  output Real hu_dyn[n](start=0.04*ones(n));

  Geometry geometry(n=n);

protected
  SI.Density rho[n] = propsLiq.rho;
  SI.SurfaceTension sigma[n] = propsLiq.sigma;
  SI.DynamicViscosity eta[n] = propsLiq.eta;
  SI.Density rho_v[n] = propsVap.eta;
equation

end BaseHoldup;
