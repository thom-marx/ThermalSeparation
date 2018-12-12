within ThermalSeparation.PhaseEquilibrium;
model ConstantK "constant value for K"
 extends BasePhaseEquilibrium;
 MediumLiquid.HenryCoefficient henryCoefficient(x_l=x_l, T=T);
 SI.Pressure p_partial[nS] "saturation pressure of each component";
 parameter Real K_user[nS];
equation
    for i in 1:nS loop
      K[i] = K_user[i];// factor_K[mapping[i,2]]* p_sat[mapping[i,2]]/p;
         p_partial[i]  =x_l[mapping[i,2]].* factor_K[mapping[i,2]].* p_sat[mapping[i,2]];
    end for;

    /*** saturation pressure of a mixture at a temperature T as sum of all partial saturation pressures ***/
      p_bubble  = sum(p_partial);
end ConstantK;
