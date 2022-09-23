within ThermalSeparation.PressureLoss;
package Reboiler
  partial model BasePressureLoss
    input SI.Pressure p_in;
    input SI.Pressure p_out;
    input Real eps_liq;
    input SI.Density rho_l;
    input SI.Density rho_v;
    input SI.Diameter d_HX;
    input SI.Length length_HX;
    input SI.Diameter d_tube;
    input Real Nw;
    input Real zeta;
    output SI.VolumeFlowRate Vdot;

    SI.Pressure deltaP = p_in-p_out;
  end BasePressureLoss;

  model TubeHX "straight-tube reboiler"
    extends BasePressureLoss;
  equation
    Vdot =  sqrt(max(1e-8,(p_in - p_out - eps_liq*d_HX*Modelica.Constants.g_n*rho_l)*2*length_HX^2*(d_HX-(Nw+1)*d_tube)^2/(zeta*Nw*(eps_liq*rho_l + (1-eps_liq)*rho_v))));
  end TubeHX;

  model Linear "linear pressure loss with nominal values"
    extends BasePressureLoss;
    parameter SI.Pressure dp_nom "nominal pressure loss";
    parameter SI.VolumeFlowRate Vdot_nom "nominal volume flow rate";
  equation
    Vdot = sign(p_in - p_out)*max(0, (abs(p_in - p_out)/(dp_nom/Vdot_nom)));
  end Linear;
end Reboiler;
