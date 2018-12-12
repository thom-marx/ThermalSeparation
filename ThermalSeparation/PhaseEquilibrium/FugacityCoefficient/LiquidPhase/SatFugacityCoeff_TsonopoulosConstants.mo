within ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.LiquidPhase;
model SatFugacityCoeff_TsonopoulosConstants
                                            // Properties of Gases & Liquids (5th edition) 4.15
replaceable package Medium = 
      ThermalSeparation.Media.WaterBasedLiquid.N2_H2O_O2 
  constrainedby ThermalSeparation.Media.BaseMediumLiquid;

parameter Integer nS=2;
parameter Integer reorgLiq[nS]={1,2};
input Integer NoOfEq[nS];
parameter SI.Temperature Tcrit[nS] = {Medium.Tcrit[reorgLiq[i]] for i in 1:nS};
parameter SI.Pressure pcrit[nS] = {Medium.pcrit[reorgLiq[i]] for i in 1:nS};
parameter Units.DipoleMoment mu[nS] = {Medium.mu[reorgLiq[i]] for i in 1:nS};

parameter Real a[nS](fixed = false);
parameter Real b[nS](fixed = false);

protected
parameter Real mu_r[nS](fixed = false);
parameter Real pcrit_atm[nS]( fixed = false);
initial equation
for m in 1:nS loop
pcrit_atm[m] = pcrit[m]*9.8692*1e-6;

mu_r[m] = 10^5*mu[m]^2 * pcrit_atm[m]/Tcrit[m]^2;

  if NoOfEq[m] == 1 then
                 a[m] = 0;
                 b[m] = 0;
                           else
    if NoOfEq[m] == 2 then
                   a[m] = -2.14 * 10^(-4)*mu_r[m] - 4.308 * 10^(-21)*mu_r[m]^8;
                   b[m] = 0;
                                                                                else
      if NoOfEq[m] == 3 then
                     a[m] = -2.188 * 10^(-4) * mu_r[m]^4 - 7.831 * 10^(-21) * mu_r[m]^8;
                     b[m] = 0;
                                                                                         else
        if NoOfEq[m] == 4 then
                       a[m] = 0.0878;
                       b[m] = 0.00908 + 0.0006957*mu_r[m];
                                      else
          if NoOfEq[m] == 5 then
                         a[m] = 0.0878;
                         b[m] = 0.0525;
                                        else

                           a[m] = -0.0109;
                           b[m] = 0;

          end if;
        end if;
      end if;
    end if;
  end if;
end for;
end SatFugacityCoeff_TsonopoulosConstants;
