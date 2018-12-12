within ThermalSeparation.PhaseEquilibrium.ActivityCoefficient;
model UNIFAC
            //** Alle Werte sowie Beispiel in Gmehling & Kolbe S. 246 ff.*** // Nur für binäre Gemische
  extends BaseActivityCoefficient;

   //müssen ins Medienmodel:
   parameter Integer SG_1=2;//Anzahl der verschiedenen Strukturgruppen einer Komponente (z.B. '1'=n-Hexan, '2'=Butanon)
   parameter Integer SG_2=3;
   parameter Integer SG_v=3;//Anzahl der insgesamt verschiedenen Strukturgruppen
   parameter Real Q_m[SG_v]={0.848,0.54,1.488};//relative van der Waalsschen Gruppenoberflächen
   parameter Real Q_1[SG_1]={0.848,0.54};
   parameter Real Q_2[SG_2]={0.848,0.54,1.488};
   parameter Real R_m[SG_v]={0.9011,0.6744,1.6724};//relative van der Waalsschen Gruppenvolumina
//    parameter Real HSG_1[SG_1]={1,0,0};// '1'=Hauptgruppe ; '0' = Untergruppe
//    parameter Real HSG_2[SG_2]={0,0,1};
   //parameter Real b_1[SG_1,SG_2]={{0,1},{
   parameter Real A_m[SG_v,SG_v]={{0,0,476.4},{0,0,476.4},{26.76,26.76,0}}; //Matrix mit Werten für die Gruppenwechselwirkungsparameter der Hauptstrukturgruppen
   parameter Real A_1[SG_1,SG_1]={{0,0},{0,0}};
   parameter Real A_2[SG_2,SG_2]={{0,0,476.4},{0,0,476.4},{26.76,26.76,0}};

//    parameter Real A_21[SG_2,SG_1]={{0,0},{0,0},{26.76,26.76}};
   parameter Real v[nS,SG_v]={{2,4,0},{1,1,1}}; // Anzahl der einzelnen Strukturgruppen innerhalb eines Moleküls

protected
   Real psi_m[SG_v,SG_v];
   Real psi_1[SG_1,SG_1];
   Real psi_2[SG_2,SG_2];
//   Real psi_21[n,SG_v,SG_v];
   Real X_m[SG_v];
   Real X_1[SG_1];
   Real X_2[SG_2];
   Real theta_m[SG_v];
   Real theta_1[SG_1];
   Real theta_2[SG_2];
   Real ln_Gamma_m[SG_v];//Gruppenaktivitätskoeffizient der Mischung
   Real ln_Gamma_r_1[SG_1];//Gruppenaktivitätskoeffizient der Reinstoffe
   Real ln_Gamma_r_2[SG_2];
   Real V[nS];
   Real F[nS];
   Real r[nS];
   Real q[nS];
   Real ln_gamma_C[nS];
   Real ln_gamma_R[nS];

equation
for i in 1:nS loop
  r[i]=sum(v[i,u]*R_m[u] for u in 1:SG_v);
  q[i]=sum(v[i,u]*Q_m[u] for u in 1:SG_v);

  V[i]=r[i]/sum(r[k]*x_l[k] for k in 1:nS);
  F[i]=q[i]/sum(q[k]*x_l[k] for k in 1:nS);

 end for;

// Mischungs-Gammas

//     for j in 1:n loop
//       for i in 1:SG_v loop
//         for k in 1:SG_v loop
//           psi_m[j,i,k]=exp(-A_m[i,k]/T[j]);
//           //psi_21[j,k,i]=exp(-A_21[k,i]/T[j]);
//         end for;
//       end for;
//     end for;

 for k in 1:SG_v loop
      ln_Gamma_m[k]=Q_m[k]*(1-Modelica.Math.log(sum(theta_m[m]*(psi_m[m,k]) for m in 1:SG_v))-sum(theta_m[m]*psi_m[k,m]/sum(theta_m[n]*psi_m[n,m] for n in 1:SG_v) 
                                                                                                  for m in 1:SG_v));

      X_m[k]=sum(v[i,k]*x_l[i] for i in 1:nS)/sum(sum(v[i,u]*x_l[i] for u in 1:SG_v) 
                                                                                    for i in 1:nS);
      theta_m[k]=Q_m[k]*X_m[k]/sum(Q_m[u]*X_m[u] for u in 1:SG_v);
     for i in 1:SG_v loop
        psi_m[k,i]=exp(-A_m[k,i]/T);
        //psi_21[j,k,i]=exp(-A_21[k,i]/T[j]);
     end for;
 end for;

// Reinstoff-Gammas

  for i in 1:SG_1 loop
    ln_Gamma_r_1[i]=Q_1[i]*(1-Modelica.Math.log(sum(theta_1[m]*(psi_1[m,i]) for m in 1:SG_1))-sum(theta_1[m]*psi_1[i,m]/sum(theta_1[n]*psi_1[n,m] for n in 1:SG_1) 
                                                                                                  for m in 1:SG_1));
    X_1[i]=v[1,i]*x_l[1]/sum(v[1,u]*x_l[1] for u in 1:SG_1);
    theta_1[i]=Q_1[i]*X_1[i]/sum(Q_1[u]*X_1[u] for u in 1:SG_1);
   for k in 1:SG_1 loop
     psi_1[i,k]=exp(-A_1[i,k]/T);
   end for;
  end for;
  for i in 1:SG_2 loop
    ln_Gamma_r_2[i]=Q_2[i]*(1-Modelica.Math.log(sum(theta_2[m]*(psi_2[m,i]) for m in 1:SG_2))-sum(theta_2[m]*psi_2[i,m]/sum(theta_2[n]*psi_2[n,m] for n in 1:SG_2) 
                                                                                                  for m in 1:SG_2));
    X_2[i]=v[2,i]*x_l[2]/sum(v[2,u]*x_l[2] for u in 1:SG_2);
    theta_2[i]=Q_2[i]*X_2[i]/sum(Q_2[u]*X_2[u] for u in 1:SG_2);
   for k in 1:SG_2 loop
     psi_2[i,k]=exp(-A_2[i,k]/T);
   end for;
  end for;

     for i in 1:nS loop
       if i==1 then
            ln_gamma_R[i]=sum(v[i,k]*(ln_Gamma_m[k]-ln_Gamma_r_1[k]) for k in 1:SG_1);
       else
            ln_gamma_R[i]=sum(v[i,k]*(ln_Gamma_m[k]-ln_Gamma_r_2[k]) for k in 1:SG_2);
       end if;

            ln_gamma_C[i]=1-V[i]+Modelica.Math.log(V[i])-5*q[i]*(1-V[i]/F[i]+Modelica.Math.log(V[i]/F[i]));
            gamma[i]=exp(ln_gamma_C[i]+ln_gamma_R[i]);
     end for;

  annotation (Documentation(info="<html>
<p><h4><font color=\"#008000\">UNIFAC</font></h4></p>
<p><br/>The UNIFAC model uses interactions between the functional groups of molecules for activity estimation for non-ideal mixtures. </p>
<p>This model splits up the activity coefficient for each species in the system into two components; a combinatorial &gamma;c and a residual component &gamma;r.</p>
<p>The combinatorial part is the same as in the UNIQUAC model. The residual part is due to the interactions of the functional groups. </p>
<p>There are some parameters which appear in both parts, like the group surface area and volume contributions R and Q(given in &QUOT;References&QUOT;) and other interaction parameters which only occur in the residual part like A. A is due to the interaction energy between groups and are listed in the following references. Note that it is not the case that Anm=Amn.</p>
<p>All theses Parameters for each group interaction have to be put in the media model. Just as the number of different groups SG_v, the number of different groups within one molecule SG_nS and the number of occurrences of the functional group on each molecule v[nS,SG_v].</p>
<p><u><font style=\"color: #008000; \">References</font></u></p>
<p>UNIFAC Structural Groups and Parameters</p>
<p>http://www.aim.env.uea.ac.uk/aim/info/UNIFACgroups.html</p>
<p>H. K. Hansen, P. Rasmussen, A. Fredenslund, M. Schiller, and J Gmehling (1991) Ind. Eng. Chem. Res. 30, 2352-2355.</p>
<p>R. Wittig, J. Lohmann, and J. Gmehling (2003) Ind. Eng. Chem. Res. 42, 183-188</p>
<p>K. Balslev and J. Abildskov (2002) Ind. Eng. Chem. Res. 41, 2047-205. </p>
</html>"));
end UNIFAC;
