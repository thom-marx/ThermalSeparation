within ThermalSeparation.Media.Correlations.DiffusionCoefficient.Liquid;
package Ion
  "package containing different models to calculate the diffusion coefficient of an ion in a liquid solution"
  partial model BaseIons
    "Base class to calculate diffusion coefficient of an ion in a liquid solution"

    parameter Integer nS=3;
    parameter Integer ic[nS]=zeros(nS)
      "ionic charge of each component (i.e.: H+ = 1, HSO3- = -1, H2O = 0";

    input SI.Temperature T;
    output SI.DiffusionCoefficient D[nS];
  equation

  end BaseIons;

  model IndefDilutedIons "indefinitely diluted ions"
    //see Horvath: Handbook of aqueous electrolyte solutions,1985, p. 289
    // values for lambda0 are tabulated for example in Horvath: Handbook of aqueous electrolyte solutions,1985, p. 262
    extends BaseIons;
      parameter Real lambda0[nS]
      "limiting equivalent conductivities of ions with respect to the solvent in (cm^2 / ohm), zeros is taken for molecules";
  equation

      for i in 1:nS loop
        if ic[i] == 0 then
          D[i] = 0;
        else
          //the factor 1e-4 is the unit conversion from cm2/s in m2/s
          D[i] = lambda0[i]* Modelica.Constants.R * T / (abs(ic[i])*Modelica.Constants.F^2)*1e-4;
        end if;
      end for;

  end IndefDilutedIons;

  model Const_Ion
    extends ThermalSeparation.Media.Correlations.DiffusionCoefficient.Liquid.Ion.BaseIons;
    parameter SI.DiffusionCoefficient D_const[nS]=fill(1e-9,nS);
  equation
    D=D_const;
  end Const_Ion;
end Ion;
