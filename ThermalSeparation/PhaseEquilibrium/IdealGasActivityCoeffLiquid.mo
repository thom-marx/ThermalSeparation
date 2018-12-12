within ThermalSeparation.PhaseEquilibrium;
model IdealGasActivityCoeffLiquid
  "gas is ideal, non-ideality of liquid taken into account using an activity coefficient model"
 extends BasePhaseEquilibrium;
 MediumLiquid.HenryCoefficient henryCoefficient(x_l=x_l, T=T);
 MediumLiquid.ActivityCoefficient activityCoeff(T=T,x_l=x_l);
 SI.Pressure p_partial[nS] "saturation pressure of each component";
equation
   for i in 1:nS loop
    if MediumLiquid.henry[mapping[i,2]] then
      K[i] = factor_K[mapping[i,2]]*henryCoefficient.He[mapping[i,2]]*activityCoeff.gamma[mapping[i,2]]/p;
        p_partial[i]  =x_l[mapping[i,2]].* factor_K[mapping[i,2]].* henryCoefficient.He[mapping[i,2]].*activityCoeff.gamma[mapping[i,2]];
    else
      K[i] = factor_K[mapping[i,2]]* p_sat[mapping[i,2]]*activityCoeff.gamma[mapping[i,2]]/p;
        p_partial[i]  =x_l[mapping[i,2]].* factor_K[mapping[i,2]].* p_sat[mapping[i,2]].*activityCoeff.gamma[mapping[i,2]];
    end if;
   end for;

   /*** saturation pressure of a mixture at a temperature T as sum of all partial saturation pressures ***/
      p_bubble  = sum(p_partial);
end IdealGasActivityCoeffLiquid;
