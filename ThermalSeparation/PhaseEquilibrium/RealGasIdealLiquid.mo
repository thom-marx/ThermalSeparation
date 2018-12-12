within ThermalSeparation.PhaseEquilibrium;
model RealGasIdealLiquid
  "gas fugacity coefficient not equal to one, liquid ideal"
 extends BasePhaseEquilibrium;
  MediumLiquid.HenryCoefficient henryCoefficient(x_l=x_l, T=T);
  MediumVapour.FugacityCoefficient fugacityCoeff(T=T, p=p, x=x_v, v=v_v);
  SI.Pressure p_partial[nS] "saturation pressure of each component";
equation
     for i in 1:nS loop
    if MediumLiquid.henry[i] then
      K[i] = factor_K[mapping[i,2]]*henryCoefficient.He[mapping[i,2]]/(p*fugacityCoeff.phi[mapping[i,1]]);
      p_partial[i]  =x_l[mapping[i,2]].* factor_K[mapping[i,2]].* henryCoefficient.He[mapping[i,2]]./fugacityCoeff.phi[mapping[i,1]];
    else
      K[i] = factor_K[mapping[i,2]]* p_sat[mapping[i,2]]/(p*fugacityCoeff.phi[mapping[i,1]]);
      p_partial[i]  =x_l[mapping[i,2]].* factor_K[mapping[i,2]].* p_sat[mapping[i,2]]./fugacityCoeff.phi[mapping[i,1]];
    end if;
    end for;

   /*** saturation pressure of a mixture at a temperature T as sum of all partial saturation pressures ***/
      p_bubble  = sum(p_partial);
end RealGasIdealLiquid;
