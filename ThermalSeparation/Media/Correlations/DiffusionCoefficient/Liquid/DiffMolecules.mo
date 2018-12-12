within ThermalSeparation.Media.Correlations.DiffusionCoefficient.Liquid;
model DiffMolecules
  "liquid binary MS diffusion coefficients in a solution which does not contain ions"
  //calculates the binary MS diffusion coefficient based on the binary infinite dilution diffusion coefficient
  //Gleichung zur Bestimmung von D auf richtige Implementierung getestet für nS=2 und nS=3 (n=1)
  parameter Integer nS(min=2)=2;
  parameter Integer ic[nS]=zeros(nS) "electric charge";
  parameter Boolean has_etaSubstance[nS];

   input SI.Temperature T;
   input SI.Pressure p;
   input SI.MoleFraction x[nS];
   input SI.DynamicViscosity eta[nS]
    "viscosity of pure liquid at system pressure and temperature";

  //The values on the diagonal do not have a physical meaning, nevertheless they are part of the equation to calculate D
  //(The expression if first multiplied and than divided by those values.)
 SI.DiffusionCoefficient D0[nS,nS]= d_molecules.D
    "matrix of binary infinite dilution diffusion coefficients";

   //Diff.koeff. z.B.: D12, D13, D14, D23, D24, D34
   output SI.DiffusionCoefficient D[a];

/*** Diffusion model for molecules (=nonelectrolytes or undissociated electrolytes) ***/
   replaceable model D_Molecules =
      ThermalSeparation.Media.Correlations.DiffusionCoefficient.Liquid.Molecule.Const_Molecule
                                               constrainedby
    ThermalSeparation.Media.Correlations.DiffusionCoefficient.Liquid.Molecule.BaseMolecules
                                                              annotation(choicesAllMatching=true);
   D_Molecules d_molecules( nS=nS, a=a, T=T, p=p, x=x, eta=eta, has_etaSubstance = has_etaSubstance);

protected
   parameter Integer counter1[nS-1]={nS-i+1 for i in 2:nS};
   Integer counter[nS];  //für nS=4: {0,3,2,1};
    SI.MoleFraction x_mod[nS]
    "modified mole fraction to treat liquid systems with dissolved ideal gas components";
    SI.MoleFraction x_mod0[nS]
    "modified mole fraction where the x of components with no liquid viscosity are set to zero";
    parameter Boolean eta_pure(fixed=false)
    "true, if for every component a liquid viscosity can be provided";

  parameter Integer aux[:] = {1,3,6,10,15, 21, 28, 36, 45};
  parameter Integer a = aux[nS-1]
    "number of binary diffusion coefficients depending on the number of substances";

equation
 counter[1]=0;
 counter[2:nS] = counter1;
 if eta_pure then
   x_mod=x;
   x_mod0 = x "x_mod0 not used in this case";
 else
   for i in 1:nS loop
     x_mod0[i] = if has_etaSubstance[i] then x[i] else 0;
     x_mod[i] = x_mod0[i]/sum(x_mod0[:]);
   end for;
 end if;

    for i in 1:nS-1 loop
      for k in i+1:nS loop
        /*** MOLECULES ***/
           //Kooijman, H. A. and Taylor, R.: Estimation of Diffusion Coefficients in Multicomponent Liquid Systems: equation (23)
          D[k-i+sum(counter[1:i])] = D0[i,k]^x_mod[k] * D0[k,i]^x_mod[i] * product((D0[i,:].*D0[k,:]).^(x_mod[:]/2)) /((D0[i,i]*D0[k,i])^(x_mod[i]/2)) /((D0[i,k]*D0[k,k])^(x_mod[k]/2));
      end for;
      end for;

initial algorithm
  eta_pure:=true;
  for i in 1:nS loop
    if not has_etaSubstance[i] then
             eta_pure:=false;
    end if;
    end for;

  annotation (Documentation(info="<html>
<p>[1] Kooijman, H. A. and Taylor, R.: Estimation of Diffusion Coefficients in Multicomponent Liquid Systems, Ind. Eng. Chem. Res. 1991, 30, pp. 1217-1222</p>
</html>"));
end DiffMolecules;
