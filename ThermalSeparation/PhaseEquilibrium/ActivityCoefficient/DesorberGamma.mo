within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient;
model DesorberGamma
  "frickel gamma des Desorbers (nur für CO2-H2O Modell) + komische Definition von xGas"
  extends BaseActivityCoefficient;
equation

//der Term (1+x_l[j,1]) beinhaltet die etwas seltsame Definition von xGas --> in gamma reingepackt
//damit man das nicht in BaseStage ändern muß (dann würde das ja für alles gelten)

       gamma[1] = 1*(exp(exp(8.99511377 - 0.0203825*T + 2.23e-5*T^2)*x_l[1]))/(1+x_l[1]);
        //  gamma_aux[j,1] = exp(exp(8.99511377 - 0.0203825*T[j] + 2.23e-5*T[j]^2)*x_l[j,1]);
    gamma[2]=1;
    gamma[3]=1;

end DesorberGamma;
