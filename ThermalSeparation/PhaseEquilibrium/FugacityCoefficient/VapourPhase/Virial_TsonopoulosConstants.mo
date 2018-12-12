within ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.VapourPhase;
model Virial_TsonopoulosConstants
                                            // Properties of Gases & Liquids 4.15
replaceable package MediumVapour = 
      ThermalSeparation.Media.IdealGasMixtures.N2_H2O_O2_CO2_SO2_HCl_HF 
                                                            constrainedby
    Media.BaseMediumVapour;

parameter Integer nS=2;
parameter Integer reorgVap[nS] = {1,2};
input Integer NoOfEq[nS];
parameter SI.Temperature Tcrit[nS] = {MediumVapour.Tcrit[reorgVap[i]] for i in 1:nS};
parameter SI.Pressure pcrit[nS] = {MediumVapour.pcrit[reorgVap[i]] for i in 1:nS};
parameter Units.DipoleMoment my[nS] = {MediumVapour.mu[reorgVap[i]] for i in 1:nS};

parameter Real a[nS](fixed = false);
parameter Real b[nS](fixed = false);

protected
parameter Real my_r[nS](fixed = false);
parameter Real pcrit_atm[nS]( fixed = false);
initial equation
for m in 1:nS loop

pcrit_atm[m] = pcrit[m]*9.8692*1e-6;

my_r[m] = 10^5*my[m]^2 * pcrit_atm[m]/Tcrit[m]^2;

  if NoOfEq[m] == 1 then
                 a[m] = 0;
                 b[m] = 0;
                           else
    if NoOfEq[m] == 2 then
                   a[m] = -2.14 * 10^(-4)*my_r[m] - 4.308 * 10^(-21)*my_r[m]^8;
                   b[m] = 0;
                                                                                else
      if NoOfEq[m] == 3 then
                     a[m] = -2.188 * 10^(-4) * my_r[m]^4 - 7.831 * 10^(-21) * my_r[m]^8;
                     b[m] = 0;
                                                                                         else
        if NoOfEq[m] == 4 then
                       a[m] = 0.0878;
                       b[m] = 0.00908 + 0.0006957*my_r[m];
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
  annotation (Documentation(info="<html>
<p><h4>Tsonopoulos Constants for Second Virial Coefficient</h4></p>
<p><u>Literature:</u> Properties of Gases &AMP; Liquids 5th Edition ( 4.13, 5.10)</p>
<p>Equation for estimation of a and b can be found in the above mentioned Literature. </p>
<p><br/>mu_r = 10^5 * mu^2 *pcrit / Tcrit^2</p>
<p><ul>
<li>mu [ debye]</li>
<li>pcrit [atm]</li>
<li>Tcrit[K]</li>
</ul></p>
</html>"));
end Virial_TsonopoulosConstants;
