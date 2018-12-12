within ThermalSeparation.HeatAndMassTransfer.ThermodynamicFactorLiq;
model UniquacBinary
//nur zur Überprüfung von UNIQUACMulticomponent

  parameter Integer n = 1;
  parameter Integer nS= 2;
  constant Real R=Modelica.Constants.R;
  input SI.Temperature T[n];
  input SI.MoleFraction x[n,nS];
  output Real Gamma[n,nS-1,nS-1];

protected
  parameter Real q[nS]=fill(1,nS);
  parameter Real r[nS]=fill(1,nS);
  parameter Real lambda[nS,nS]=fill(1,nS,nS);
  parameter Integer z=10;
  Real r2[n];
  Real q2[n];
  Real Phi[n,nS];
  Real Qc11[n];
  Real Qr11[n];
  Real Qc12[n];
  Real Qr12[n];
  Real tau12[n];
  Real tau21[n];
  Real S[n,nS];
equation
  for j in 1:n loop
  r2[j]=x[j,1]*r[1]+x[j,2]*r[2];
  q2[j]=x[j,1]*q[1]+x[j,2]*q[2];
  Phi[j,1]=x[j,1]*q[1]/q2[j];
  Phi[j,2]=x[j,2]*q[2]/q2[j];
  tau12[j]=exp(-(lambda[1,2]-lambda[1,1])/R*T[j]);
  tau21[j]=exp(-(lambda[2,1]-lambda[2,2])/R*T[j]);
  S[j,1]=Phi[j,1]+Phi[j,2]*tau12[j];
  S[j,2]=Phi[j,2]+Phi[j,1]*tau21[j];
  Qc11[j]=-2*r[1]/r2[j]+(r[1]/r2[j])^2*(x[j,1]+x[j,2])-z/2*q2[j]*(r[1]/r2[j]-q[1]/q2[j])^2;
  Qr11[j]=q[1]^2*(1-2/S[j,1]+Phi[j,1]/S[j,1]^2+Phi[j,2]*tau12[j]^2/S[j,2]^2)/q2[j];
  Qc12[j]=-r[1]/r2[j]-r[2]/r2[j]+(r[1]*r[2]/r2[j]^2)*(x[j,1]+x[j,2])-z/2*q2[j]*(r[2]/r2[j]-q[2]/q2[j])*(r[1]/r2[j]-q[1]/q2[j]);
  Qr12[j]=q[1]*q[2]*(1-tau12[j]/S[j,2]-tau21[j]/S[j,1]+Phi[j,1]*tau21[j]/S[j,1]^2+Phi[j,2]*tau12[j]/S[j,2]^2)/q2[j];
  Gamma[j,1,1]=1+x[j,1]*(Qc11[j]+Qr11[j]-Qc12[j]-Qr12[j]);
  end for;

end UniquacBinary;
