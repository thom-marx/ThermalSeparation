﻿within ThermalSeparation.FilmModel.RandomPackedColumn;
model MS "Maxwell-Stefan mass transfer - no film reaction"
extends BaseNonEqType;
extends ThermalSeparation.FilmModel.BaseClasses.MaxwellStefan(redeclare record BaseGeometry=Geometry);

 final parameter Integer aux[  :] = {1,3,6,10,15, 21, 28, 36, 45};
Geometry geometry(n=n);

      replaceable model InterfacialArea =
      ThermalSeparation.HeatAndMassTransfer.InterfacialArea.RandomPackedColumn.Constant                 constrainedby ThermalSeparation.HeatAndMassTransfer.InterfacialArea.RandomPackedColumn.BasePacked
                                                                                                      annotation(choicesAllMatching=true,Dialog( tab="Heat and Mass Transfer", group="Contact area between liquid and vapour phase"));
    InterfacialArea interfacialArea(
      redeclare record Geometry=Geometry,
      final n=n,rho=propsLiq.rho, sigma=propsLiq.sigma, eta=propsLiq.eta, Vdot_l=Vdot_l, eps_liq=eps_liq);

replaceable model MassTransferCoeffVap =
      ThermalSeparation.HeatAndMassTransfer.RandomPackedColumn.Vapour.Constant
   constrainedby ThermalSeparation.HeatAndMassTransfer.RandomPackedColumn.Vapour.BaseVapMT
     annotation(choicesAllMatching=true,Dialog( tab="Heat and Mass Transfer", group="Liquid and vapour mass transfer"));

 MassTransferCoeffVap massTransferCoeffVap(n=n, n_k=aux[nSV-1],D=D_vap, eta=propsVap.eta, sigma=propsLiq.sigma, rho=propsVap.rho, w_sup_l=w_sup_l, w_sup_v=w_sup_v, eps_liq=eps_liq,
 redeclare record Geometry =  Geometry);

replaceable model MassTransferCoeffLiq =
      ThermalSeparation.HeatAndMassTransfer.RandomPackedColumn.Liquid.Constant
   constrainedby ThermalSeparation.HeatAndMassTransfer.RandomPackedColumn.Liquid.BaseLiqMT
     annotation(choicesAllMatching=true,Dialog( tab="Heat and Mass Transfer", group="Liquid and vapour mass transfer"));

 MassTransferCoeffLiq massTransferCoeffLiq(n=n, n_k=aux[nSL-1],D=D_liq, eta=propsLiq.eta, sigma=propsLiq.sigma, rho=propsLiq.rho, w_sup=w_sup_l, eps_liq=eps_liq,
 redeclare record Geometry =  Geometry);

   replaceable model HeatTransferVapour =
      ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient.Vapour.Constant
     constrainedby ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient.Vapour.BaseVapour
                                                                                                        annotation(choicesAllMatching=true,Dialog( tab="Heat and Mass Transfer", group="Heat transfer between liquid and vapour phase"));

   HeatTransferVapour heatTransferCoeffVap(
   redeclare replaceable package Medium =  MediumVapour,
   n=n, nS=nSV,  k_av=k_av_vap, D_av=D_av_vap,
   rho=propsVap.rho,  c=c_v, props=propsVap);

     replaceable model HeatTransferLiquid =
      ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient.Liquid.Constant
     constrainedby ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient.Liquid.BaseLiquid
                                                                                                        annotation(choicesAllMatching=true,Dialog( tab="Heat and Mass Transfer", group="Heat transfer between liquid and vapour phase"));

   HeatTransferLiquid heatTransferCoeffLiq(
    redeclare replaceable package Medium = MediumLiquid,
    n=n, nS=nSL, k_av=k_av_liq, D_av=D_av_liq,
    rho=propsLiq.rho,   c=c_l, props=propsLiq);
ThermalSeparation.Units.CoefficentOfMassTransfer k_av_liq[n];
ThermalSeparation.Units.CoefficentOfMassTransfer k_av_vap[n];
SI.DiffusionCoefficient D_av_liq[n];
SI.DiffusionCoefficient D_av_vap[n];

       MediumLiquid.DiffusionCoefficient[n] diffCoeffLiq(
    T=T_l,
    p=p_v[1:n],
    x=x_l,
    eta=eta_comp);
     SI.DiffusionCoefficient D_liq[n,aux[nSL-1]]=diffCoeffLiq.D;
       MediumVapour.DiffusionCoefficient[n] diffCoeffVap(T=T_v, p=p_v[1:n]);
     SI.DiffusionCoefficient D_vap[n,aux[nSV-1]]=diffCoeffVap.D;
     SI.Velocity w_sup_l[n] "superficial liquid velocity";
     SI.Velocity w_sup_v[n] "superficial vapour velocity";
equation

           //Erstellen der K Matrix:
       for j in 1:n loop
    //Erstellen einer Matrix aus den binären Stoffübergangskoeffizienten
     for i in 1:nSL loop
       for m in i:nSL loop
         if i==m then
           //die Einträge auf der Diagonalen, werden auf irgendeinen Wert gesetzt, da eh nie benutzt
           k_l[j,i,m] = 10;
         else
           //die Einträge des Vektors werden an die richtigen Stellen auf der Matrix verteilt
           k_l[j,i,m] = max(massTransferCoeffLiq.k[j,m-2+i],1e-10);
           end if;
         end for;
         for m in 1:i-1 loop
           //k_12 = k_21, k_13 = k_31 etc.
           k_l[j,i,m] = max(k_l[j,m,i],1e-10);
         end for;
     end for;
         end for;
                for j in 1:n loop
    //Erstellen einer Matrix aus den binären Stoffübergangskoeffizienten
     for i in 1:nSV loop
       for m in i:nSV loop
         if i==m then
           //die Einträge auf der Diagonalen, werden auf irgendeinen Wert gesetzt, da eh nie benutzt
           k_v[j,i,m] = 10;
         else
           //die Einträge des Vektors werden an die richtigen Stellen auf der Matrix verteilt
           k_v[j,i,m] = max(massTransferCoeffVap.k[j,m-2+i],1e-10);
           end if;
         end for;
         for m in 1:i-1 loop
           //k_12 = k_21, k_13 = k_31 etc.
           k_v[j,i,m] = max(k_v[j,m,i],1e-10);
         end for;
     end for;
      end for;

   alphaLiq=heatTransferCoeffLiq.alpha;
   alphaVap=heatTransferCoeffVap.alpha;

 for j in 1:n loop
     D_av_vap[j]=sum(D_vap[j,:])/aux[nSV-1];
      D_av_liq[j]=sum(D_liq[j,:])/aux[nSL-1];
    k_av_vap[j]=sum(massTransferCoeffVap.k[j,:])/aux[nSV-1];
      k_av_liq[j]=sum(massTransferCoeffLiq.k[j,:])/aux[nSL-1];
     A_I[j] = interfacialArea.a[j]*geometry.H*geometry.A/n;
  end for;

     /*** superficial velocity ***/
         for j in 1:n loop
  w_sup_l[j] = max(1e-10,Vdot_l[j]/geometry.A);
    w_sup_v[j] = max(1e-10,Vdot_v[j]/geometry.A);
      end for;
end MS;
