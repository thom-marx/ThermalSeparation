within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient;
model DesorberGammaBremsen
  "frickel gamma des Desorbers (nur für CO2-H2O Modell) + komische Definition von xGas + bremsen"
  extends BaseActivityCoefficient;
Real bremsen(start=0);
Real K=1;
Real bremsenStart=0.001;
parameter Real bremsen_start= 0;
initial equation
  bremsen=bremsen_start;
equation

//der Term (1+x_l[j,1]) beinhaltet die etwas seltsame Definition von xGas --> in gamma reingepackt
//damit man das nicht in BaseStage ändern muß (dann würde das ja für alles gelten)
//bremsen[j] = (max(1/0.0001,1/x_l_total[j,3])^0.3)-10000^0.3 ;
der(bremsen)=K* ((max(1/bremsenStart,1/x_l[3])^1)-(1/bremsenStart)^1 - bremsen);
       gamma[1] = 1*(exp(exp(8.99511377 - 0.0203825*T + 2.23e-5*T^2)*x_l[1]))/(1+x_l[1])+bremsen;
        //  gamma_aux[j,1] = exp(exp(8.99511377 - 0.0203825*T[j] + 2.23e-5*T[j]^2)*x_l[j,1]);
    gamma[2]=1;
    gamma[3]=1;

end DesorberGammaBremsen;
