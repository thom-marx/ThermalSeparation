within ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.VapourPhase;
model Virial
/*** see for example Reid, Prausnitz: The Properties of Gases & Liquids, 4th edition, p. 145 and p. 79 ***/
 extends BaseFugacityCoefficient;

/*** second virial coefficients for mixtures ***/

parameter Integer NoOfEq[nS] = {MediumVapour.eq_Tsonopoulos[reorgVap[i]] for i in 1:nS};

ThermalSeparation.PhaseEquilibrium.FugacityCoefficient.VapourPhase.Virial_SecondVirialCoefficient
    svc(
    n=n,
    nS=nS,
    reorgVap=reorgVap,
    redeclare replaceable package MediumVapour = MediumVapour,
    T=T,
    NoOfEq=NoOfEq);

protected
  Real Bij[n,nS,nS]=svc.B;
  Real B[n];
  Real Bi[n,nS];
  Real c[n,nS] "Hilfsvariable ohne physik. Bedeutung";

equation
/*** determination of the second virial coefficient ***/

for m in 1:n loop
  B[m]=sum(sum(y[m,i]*y[m,j]*Bij[m,i,j] for i in 1:nS) for j in 1:nS);
  for i in 1:nS loop
    c[m,i]=sum(y[m,j]*Bij[m,i,j] for j in 1:nS);
    Bi[m,i]= 2*c[m,i]-B[m];
    phi_aux[m,i]=exp((p[m]*Bi[m,i])/(Modelica.Constants.R*T[m]));
  end for;
end for;

    annotation (Documentation(info="<html>
<p align=\"justify\"><h4>Virial Equation of State </h4></p>
<p>phi_vap ist durch die Gleichung</p>
<p><br/>phi_vap_i = exp [ 2*(sum(yj*Bij) - B) * P/ (R*T)]</p>
<p>gegeben. [2]</p>
<p>Bij sind die Virialkoeffizienten der Reinstoffe bzw. die sogenannten Cross-Koeffizienten.</p>
<p>Bei der Berechnung der Cross-Koeffizienten treten bin&auml;re Interaktionskoeffizienten auf. Diese sind oft nicht viel gr&ouml;&szlig;er/kleiner als 0.</p>
<p>Einige k&ouml;nnen in &QUOT;Fluid Phase Equilibria 260 (2007) 354-358&QUOT; nachgelesen werden.</p>
<p>Desweiteren werden f&uuml;r die Berechnung der Cross-Koeffizienten Mischungsgr&ouml;&szlig;en f&uuml;r Tcrit, pcrit, omega und vcrit ben&ouml;tigt.</p>
<p>Diese findet man in Lit [1, 5.10] oder in Lit[4].</p>
<p>Die Berechnung erfolgt mit Hilfe der Tsonopoulos Konstanten. Diese h&auml;ngen vom Dipolmoment ab. [1]</p>
<p>B ist der zweite Koeffizient der Virialgleichung und ergibt sich aus:</p>
<p><br/>B = sum(sum(yi*yj*Bij))</p>
<p><br/><br/>References:</p>
<p>[1] Poling, Prausnitz, O&apos;Connell : The Properties of Gases &AMP; Liquids 5th Edition [4.13, 5.10]</p>
<p>[2] Gmehling, Kolbe: Thermodynamik [p. 120]</p>
<p>[3] Dymond, Smith : The Virial Coefficients of Pure Gases and Mixtures</p>
<p>[4] Fluid Phase Equilibria 260 (2007) 354-358</p>
</html>"));
end Virial;
