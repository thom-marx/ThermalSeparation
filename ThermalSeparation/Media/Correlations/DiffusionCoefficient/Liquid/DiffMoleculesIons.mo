within ThermalSeparation.Media.Correlations.DiffusionCoefficient.Liquid;
model DiffMoleculesIons
  "liquid binary MS diffusion coefficients for a solution which may contain ions"
  //calculates the binary MS diffusion coefficient based on the binary infinite dilution diffusion coefficient
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

   /*** Diffusion model for ions ***/
   replaceable model D_Ions =
      ThermalSeparation.Media.Correlations.DiffusionCoefficient.Liquid.Ion.Const_Ion
    constrainedby
    ThermalSeparation.Media.Correlations.DiffusionCoefficient.Liquid.Ion.BaseIons
                                                                 annotation(choicesAllMatching=true);
   D_Ions d_ion( nS=nS, T=T, ic=ic);

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

        /*** IONS ***/
         if ic[i] <> 0 and ic[k] <>0 then
         //The binary diffusion coefficient of an ion i with respect to an ion j is set to the mean of the effective diffusivities of the two ions:
          D[k-i+sum(counter[1:i])] = (d_ion.D[i] + d_ion.D[k])/2;
           elseif ic[i] <> 0 and ic[k] == 0 then
          //The binary diffusion coefficient of the ion with respect to any molecular species is set equal to the diffusivity of the ion in the liquid mixture:
          D[k-i+sum(counter[1:i])] = d_ion.D[i];
         elseif ic[i] ==0 and ic[k] <> 0 then
         D[k-i+sum(counter[1:i])] = d_ion.D[k];
         else

        /*** MOLECULES ***/
           //Kooijman, H. A. and Taylor, R.: Estimation of Diffusion Coefficients in Multicomponent Liquid Systems: equation (23)
          D[k-i+sum(counter[1:i])] = D0[i,k]^x_mod[k] * D0[k,i]^x_mod[i] * product((D0[i,:].*D0[k,:]).^(x_mod[:]/2)) /((D0[i,i]*D0[k,i])^(x_mod[i]/2)) /((D0[i,k]*D0[k,k])^(x_mod[k]/2));
        end if;
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
end DiffMoleculesIons;
