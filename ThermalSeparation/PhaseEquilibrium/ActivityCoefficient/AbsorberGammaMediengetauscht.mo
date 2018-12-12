within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient;
model AbsorberGammaMediengetauscht
  "frickel gamma des Absorbers (nur für CO2-H2O Modell) + komische Definition von xGas, getauschte Medien"
  extends BaseActivityCoefficient;

equation
//der Term (1+x_l[j,1]) beinhaltet die etwas seltsame Definition von xGas --> in gamma reingepackt
//damit man das nicht in BaseStage ändern muß (dann würde das ja für alles gelten)
       gamma[3] = (exp(exp(8.99511377 - 0.0203825*T + 2.23e-5*T^2)*x_l[2]))/(1+x_l[2]);
        //  gamma_aux[j,1] = exp(exp(8.99511377 - 0.0203825*T[j] + 2.23e-5*T[j]^2)*x_l[j,1]);
    gamma[1]=1;
    gamma[2]=1;

end AbsorberGammaMediengetauscht;
