within ThermalSeparation.HeatAndMassTransfer.DiffusionCoeff.Liquid;
model test
  //Calculates the matrix of multicomponent mass transfer coefficients [k] using the Fick diffusion matrix [D],
  //method of Toor, Stewart and Prober (Taylor: 8.8.2, p. 214)
  //Validation using Example 8.8.1
  parameter Integer n= 3;
  parameter Real[:,:] D=[2.10174, -0.117519; -0.089253, 2.11751]*1e-6;
   //parameter Real[:,:] D=[2.093, -0.13; -0.067,  2.148]*1e-6;
  Real eigD[n-1] = Modelica.Math.Matrices.LAPACK.dgeev_eigenValues(D);
  Real eigSh[n-1];
  Real eigk[n-1];
  Real eigSc[n-1];
  parameter Real[:,:] I=diagonal(ones(n-1));

 Real prod[n-1];
 Real prod2[n-1,n-1,n-1];

     Integer m[n-1];
      Real test[n-1,n-2];
        Integer m2[n-1];
  Real D2[  n-1,n-2,n-1,n-1];

  Real k_help[n-1,n-1,n-1];
  Real k[n-1,n-1];

equation
  for i in 1:n-1 loop
    eigSc[i]=8.819e-6/2.81/eigD[i];
    eigSh[i] = 0.023*21752^0.83 *(eigSc[i])^0.44;
    eigk[i] = eigSh[i]*eigD[i]/2.21e-2;
  end for;

  //Problem: D ist nicht konstant sondern abh‰ngig von T und c, deswegen kann es sein, daﬂ man hier gar keinen algorithm verwenden darf!
algorithm
  for i in 1:n-1 loop
    m[i]:=0;
    for j in 1:n-1 loop
      if i==j then
        m[i]:=1;
      else
        test[i,j-m[i]]:=eigD[i] - eigD[j];
      end if;
      end for;
      prod[i] :=product(test[i, :]);
      end for;

        for i in 1:n-1 loop
        m2[i]:=0;
    for j in 1:n-1 loop
      if i==j then
        m2[i]:=1;
        else
    D2[i,j-m2[i],:,:]:=D - eigD[j]*I;
    end if;
      end for;

      end for;

      for i in 1:n-1 loop
        prod2[i,:,:]:=I;
        for j in 1:n-2 loop
      prod2[i,:,:]:=prod2[i, :, :]*D2[i, j, :, :];
      end for;
        end for;

        for i in 1:n-1 loop
        k_help[i,:,:]:=eigk[i]*prod2[i, :, :]/prod[i];
        end for;

        for i in 1:n-1 loop
          k[:,:] :=k[:, :] + k_help[i, :, :];
          end for;

end test;
