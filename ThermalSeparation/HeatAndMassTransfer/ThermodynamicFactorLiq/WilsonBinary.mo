within ThermalSeparation.HeatAndMassTransfer.ThermodynamicFactorLiq;
model WilsonBinary
  //inputs:
  parameter Real T=400;
  parameter Real R=10;
  parameter Real x[2]={0.5,0.5};
  parameter Units.DipoleMoment mu[2,2]={{3,6},{8,2}};
  output Real Gamma;

protected
  parameter Real V[2]={100,10};
  Real Q[2];
  Real S[2];
  Real lambda[2];
  Real RT=R*T;
equation
  lambda[1]=(V[2]/V[1])*exp(-(mu[1,2]-mu[1,1])/RT);
  lambda[2]=(V[1]/V[2])*exp(-(mu[2,1]-mu[2,2])/RT);
  S[1]=x[1]+x[2]*lambda[1];
  S[2]=x[2]+x[1]*lambda[2];
  Q[1]=-2/S[1]+x[1]/(S[1]^2)+x[2]*lambda[2]^2/S[2]^2;
  Q[2]=-lambda[1]/S[1]-lambda[2]/S[2]+x[1]*lambda[1]/S[1]^2+x[2]*lambda[2]/S[2]^2;
  gamma=1+x[1]*(Q[1]-Q[2]);
end WilsonBinary;
