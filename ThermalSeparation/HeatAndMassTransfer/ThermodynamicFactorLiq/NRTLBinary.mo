within ThermalSeparation.HeatAndMassTransfer.ThermodynamicFactorLiq;
model NRTLBinary
  //nur zur Überprüfung von NRTLMulticomponent
  parameter Integer n=3;
  parameter Real x[3,2]={{0.8,0.2},{0.2,0.8},{0.5,0.5}};
  output Real Gamma[n];
  parameter Real alpha=2;
  parameter Real g[2,2]={{1,2},{3,4}};
  parameter Real R=10;
  parameter Real T[n]={400,300,200};
protected
  Real tau[n,2];
  Real G[n,2];
  Real S[n,2];
  Real a[n];
  Real b[n];
equation
  for j in 1:n loop
    tau[j,1]=(g[1,2]-g[1,1])/(R*T[j]);
    tau[j,2]=(g[2,1]-g[2,2])/(R*T[j]);
    G[j,1]=exp(-alpha*tau[j,1]);
    G[j,2]=exp(-alpha*tau[j,2]);

    S[j,1]=x[j,1]*x[j,2]*G[j,2];
    S[j,2]=x[j,2]*x[j,1]*G[j,1];
    a[j]=1/S[j,2]^3;
    b[j]=1/S[j,1]^3;
    Gamma[j]=1-2*x[j,1]*x[j,2]*(tau[j,2]*G[j,2]^2/(a[j])+tau[j,1]*G[j,1]^2/b[j]);
  end for;

end NRTLBinary;
