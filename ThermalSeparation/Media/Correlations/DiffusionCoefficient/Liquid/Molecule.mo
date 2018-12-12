within ThermalSeparation.Media.Correlations.DiffusionCoefficient.Liquid;
package Molecule
  "Binary diffusion coefficients for molecules or non-dissociated salts"
  partial model BaseMolecules

    parameter Integer nS=3;
    parameter Integer a=3;
      parameter Boolean has_etaSubstance[nS];

     input SI.Temperature T;
     input SI.Pressure p;
     input SI.MoleFraction x[nS];
     input SI.DynamicViscosity eta[nS];
    //input SI.Density rho_l[n];
     output SI.DiffusionCoefficient D[nS,nS];
  equation

  end BaseMolecules;

  model Const_Molecule
    "constant and equal diffusion coefficients for molecules"
    extends
      ThermalSeparation.Media.Correlations.DiffusionCoefficient.Liquid.Molecule.BaseMolecules;
    parameter SI.DiffusionCoefficient D_const_liquid = 2e-9;
  equation
    D=fill(D_const_liquid,nS,nS);
  end Const_Molecule;

  model Stokes_Einstein "Stokes-Einstein "
    extends BaseMolecules;

    parameter Real r[nS] "radius of the diffusing molecule"; //order of magnitude: =4e-13

  equation
  for i in 1:nS loop
    for k in 1:nS loop
       if i==k then
          D[i,k]=1e-9 "value has no effect on result, see class DiffCoeffLiq";
        else
          if has_etaSubstance[k] then
            //factor 1000 is necessary, since the unit of eta has to be converted
            D[i,k] = Modelica.Constants.k * T /(6*Modelica.Constants.pi*eta[k]*1000*r[i]);
          else
            //If the component which is to be the solvent does not provide a liquid viscosity, the Stoke's equation
            //cannot be used. The diffusion coefficient is set to a constant value which is however not used (see
            //class DiffCoeffLiq).
         D[i,k] = 1e-9;
          end if;
          end if;
    end for;
    end for;

    annotation (Documentation(info="<html>
<p>This correlation is only valid if the molecules of the diffusing species are very large compared to the solvent molecules.</p>
<p>Literature:</p>
<p>Taylor, R. and Krishna, R.: Multicomponent mass transfer, p. 73</p>
</html>"));
  end Stokes_Einstein;

  model Wilke_Chang_aq "Wilke-Chang for aqueous solutions"
        //Reid, Prausnitz: The Properties of Gases and Liquids, 4th edition, p. 598
    extends BaseMolecules;
    parameter Real phi[nS]
      "association factor of each substance, if this substance is to be the solvent";
    parameter SI.MolarMass MMX[nS];

  equation
  for i in 1:nS loop
    for k in 1:nS loop

        if i==k then
          D[i,k]=1e-9 "value has no effect on result, see class DiffCoeffLiq";
        else
          if has_etaSubstance[k] then
            //several factors necessary for unit conversion
            D[i,k] =7.4e-8 * (phi[k]*MMX[k]*1000)^0.5* T/(eta[k]*1000 *(MMX[i]/1000*1e6)^0.6) *1e-4;
          else
            //If the component which is to be the solvent does not provide a liquid viscosity, the Stoke's equation
            //cannot be used. The diffusion coefficient is set to a constant value which is however not used (see
            //class DiffCoeffLiq).
         D[i,k] = 1e-9;
          end if;
          end if;

    end for;
    end for;

    annotation (Documentation(info="<html>
<p>This correlation is only valid if the molecules of the diffusing species are very large compared to the solvent molecules.</p>
<p>Literature:</p>
<p>Taylor, R. and Krishna, R.: Multicomponent mass transfer, p. 73</p>
</html>"));
  end Wilke_Chang_aq;
end Molecule;
