within ThermalSeparation.Media.Correlations.DiffusionCoefficient.Vapour;
model KineticTheory
  "theoretical equation for the mutual diffusion coefficient in a low-pressure binary gas mixture (ideal gas!)"
  //z.B. Taylor, p. 68
  extends
    ThermalSeparation.Media.Correlations.DiffusionCoefficient.Vapour.BaseDiffusionCoeffGas;
    parameter Real sigma[nS];
    parameter Real epsilon_k[nS];
protected
  parameter Integer counter1[nS-1]={nS-i+1 for i in 2:nS};//für nS=4: {0,3,2,1};
  Integer counter[nS];
  Real sigma_mix[a] "Lennard-Jones force constant), constant value";
  Real epsilon_k_mix[a]
    "Lennard-Jones force constant (epsilon/k), constant value";
  final parameter Real C= 1.883e-2;

  Real omega[a] "collision integral";
  Modelica.Blocks.Tables.CombiTable1D table_omega[a](each table=collisionIntegral);

equation
   for i in 1:a loop
    table_omega[i].u = {T/epsilon_k_mix[i]};
    omega[i] = table_omega[i].y[1];
   end for;

  counter[1]=0;
  counter[2:nS] = counter1;

for i in 1:nS-1 loop
  for k in i+1:nS loop
    //the factor 1000 is due to unit conversion
 D[k-i+ sum(counter[1:i])] = ((MMX[i] + MMX[k])/(MMX[i] * MMX[k] *1000))^0.5 * C * T^(3/2)/(p*sigma_mix[k-i+ sum(counter[1:i])]^2 * omega[k-i+ sum(counter[1:i])]);
  end for;
  end for;

 for i in 1:nS-1 loop
  for k in i+1:nS loop
        sigma_mix[k-i+ sum(counter[1:i])] = 0.5*(sigma[i] + sigma[k]);
    epsilon_k_mix[k-i+ sum(counter[1:i])]=sqrt(epsilon_k[i]*epsilon_k[k]);
      end for;
  end for;

  for i in 1:nS loop
    for k in 1:nS loop
      D_matrix[i,k] = ((MMX[i] + MMX[k])/(MMX[i] * MMX[k] *1000))^0.5 * C * T^(3/2)/(p*sigma_mix[k-i+ sum(counter[1:i])]^2 * omega[k-i+ sum(counter[1:i])]);
    end for;
  end for;

  annotation (Documentation(info="<html>
<p>The diffusion coefficient is calculated by a correlation which is based on the kinetic theory of gases. The equation can be found for instance in [1] and [2]. Values for sigma and epsilon/k for the binary pairs are calculated using the values for sigma and epsilon for the pure substances. These values must be supplied in the medium model. The values for many different substances can be found in [3]. If the values for the pure substances are not available, they can be estimated using the critical data (see [2]).</p>
<p><br/>References:</p>
<p>[1] Taylor, R. and Ross, R.: Multicomponent Mass Transfer, Wiley, 1993, p. 68 </p>
<p>[2] Sherwood, T. et al.: Mass Transfer, McGraw-Hill, 1952, p. 19-23</p>
<p>[3] Svehla, R. A.: Estimated Viscosities and ..., NASA Tech. Rep. R-132, 1962</p>
</html>"));
end KineticTheory;
