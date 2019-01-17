within ThermalSeparation.Holdup;
package SprayColumn
  partial model BaseHoldup
    parameter Integer n(min=1) annotation(Dialog(enable=false));
    parameter Boolean constProp = false
      "constant values for density and surface tension";
    parameter SI.Density rho_l_const = 1000 "constant value for liquid density" annotation(Dialog(enable = constProp));
    parameter SI.Density rho_v_const = 2 "constant value for vapour density" annotation(Dialog(enable = constProp));
    parameter SI.DynamicViscosity eta_v_const = 0.075
      "constant value for surface tension"                                               annotation(Dialog(enable = constProp));
    input SI.Density rho_l[n];
    input SI.Density rho_v[n];
    input SI.DynamicViscosity eta_v[n];
    input SI.VolumeFlowRate Vdot_v[n];
    input Real eps_liq[n];
    input SI.VolumeFlowRate Vdot_l_in;
    output SI.VolumeFlowRate Vdot[n];
    output SI.Diameter d_drop[n](start=1e-3*ones(n));
    output Real n_drop[n](start=1*ones(n));

    replaceable record Geometry =
        ThermalSeparation.Geometry.SprayColumn.Geometry;
    Geometry geometry(n=n);

    SI.Density rho_liq[n] = if constProp then fill(rho_l_const,n) else rho_l;
    SI.Density rho_vap[n] = if constProp then fill(rho_v_const,n) else rho_v;
    SI.DynamicViscosity eta_vap[n] = if constProp then fill(eta_v_const,n) else eta_v;
  equation

  end BaseHoldup;

  model IdealDroplets
    extends BaseHoldup;
    parameter SI.Diameter d_drop_in=0.5e-3 "droplet diameter at liquid inlet";
  protected
    SI.Velocity c_rel[n](start=0*ones(n));
    SI.Velocity c_drop[n](start=ones(n));
    SI.Velocity c_vapour[n];
    Real zeta[n](start=ones(n));
    Real Re[n](start=100*ones(n));

  equation
    for j in 1:n loop

   // c_rel[j] =  min(200,sqrt(max(1e-5,(1-rho_v[j]/rho_l[j])*Modelica.Constants.g_n*4/3 / zeta[j] * rho_l[j]/rho_v[j] * d_drop[j])));
  c_rel[j] = ((1-rho_vap[j]/rho_liq[j])*Modelica.Constants.g_n*4/3 / 24 / max(1e-8,eta_vap[j])* rho_liq[j])^3*(geometry.A*eps_liq[j]/max(2e-5,Vdot_l_in) *d_drop_in^3)^2;
  /*** für Gleichstrom: ***/
   // c_drop[j] = max(0,abs(Vdot_v[j])/geometry.A + c_rel[j]);
  /*** für Gegenstrom: ***/
  c_drop[j] = max(0, c_rel[j] - abs(Vdot_v[j]/geometry.A));
    c_vapour[j] = abs(Vdot_v[j])/geometry.A;
    d_drop[j] = (geometry.A * geometry.H/n * eps_liq[j] / max(1e-5,n_drop[j]) *6 /Modelica.Constants.pi)^(1/3);
    Re[j] =  max(1e-5,c_rel[j] * d_drop[j]/eta_vap[j]*rho_vap[j]);
    zeta[j] = 24/Re[j];
  //  zeta[j] = 24/Re[j] + 3.73/Re[j]^0.5 - 4.82e-3*Re[j]^0.5/(1 + 3e-6*Re[j]^1.5) + 0.49;
     n_drop[j] = Vdot_l_in *6/Modelica.Constants.pi/d_drop_in^3*geometry.H/n/max(1e-5,c_drop[1]);

     Vdot[j] = c_drop[j]*geometry.A * eps_liq[j];

  end for;

  end IdealDroplets;
end SprayColumn;
