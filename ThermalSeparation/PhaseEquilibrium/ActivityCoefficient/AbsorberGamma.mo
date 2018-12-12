within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient;
model AbsorberGamma
  "frickel gamma des Absorbers (nur für CO2-H2O Modell) + komische Definition von xGas"
  extends BaseActivityCoefficient;
  parameter Real gamma_water = 1.0;
  parameter Real factor_gamma_CO2 = 1.0;
equation

//der Term (1+x_l[j,1]) beinhaltet die etwas seltsame Definition von xGas --> in gamma reingepackt
//damit man das nicht in BaseStage ändern muß (dann würde das ja für alles gelten)
       gamma[1] = factor_gamma_CO2*(exp(exp(8.99511377 - 0.0203825*T + 2.23e-5*T^2)*x_l[1]))/(1+x_l[1]);
        //  gamma_aux[j,1] = exp(exp(8.99511377 - 0.0203825*T[j] + 2.23e-5*T[j]^2)*x_l[j,1]);
    gamma[2]=gamma_water;
    gamma[3] = 1;

end AbsorberGamma;
