within ThermalSeparation.FilmModel.StructuredPackedColumn;
model EffectiveDiff
  "effective diffusion coefficient - only for diluted systems"
  import ThermalSeparation;
extends BaseNonEqType;
extends ThermalSeparation.FilmModel.BaseClasses.EffectiveDiff(redeclare record BaseGeometry =
        Geometry);

 final parameter Integer aux[  :] = {1,3,6,10,15, 21, 28, 36, 45};
Geometry geometry(n=n);

      replaceable model InterfacialArea =
      ThermalSeparation.HeatAndMassTransfer.InterfacialArea.StructuredPackedColumn.RochaStructured      constrainedby ThermalSeparation.HeatAndMassTransfer.InterfacialArea.StructuredPackedColumn.BasePacked
                                                                                                      annotation(choicesAllMatching=true,Dialog( tab="Heat and Mass Transfer", group="Contact area between liquid and vapour phase"));
    InterfacialArea interfacialArea(
      redeclare record Geometry=Geometry,
      final n=n,rho=propsLiq.rho, sigma=propsLiq.sigma, eta=propsLiq.eta, Vdot_l=Vdot_l, eps_liq=eps_liq);

replaceable model MassTransferCoeffVap =
      ThermalSeparation.HeatAndMassTransfer.StructuredPackedColumn.Vapour.RochaStructured
   constrainedby ThermalSeparation.HeatAndMassTransfer.StructuredPackedColumn.Vapour.BaseVapMT
     annotation(choicesAllMatching=true,Dialog( tab="Heat and Mass Transfer", group="Liquid and vapour mass transfer"));

 MassTransferCoeffVap massTransferCoeffVap(
 n=n, n_k=nSV,D=D_vap, eta=propsVap.eta, sigma=propsLiq.sigma, rho=propsVap.rho, w_sup_l=w_sup_l, w_sup_v=w_sup_v, eps_liq=eps_liq,
 redeclare record Geometry =  Geometry);

replaceable model MassTransferCoeffLiq =
      ThermalSeparation.HeatAndMassTransfer.StructuredPackedColumn.Liquid.RochaStructured
   constrainedby ThermalSeparation.HeatAndMassTransfer.StructuredPackedColumn.Liquid.BaseLiqMT
     annotation(choicesAllMatching=true,Dialog( tab="Heat and Mass Transfer", group="Liquid and vapour mass transfer"));

 MassTransferCoeffLiq massTransferCoeffLiq(
  n=n, n_k=nSL,D=D_liq, eta=propsLiq.eta, sigma=propsLiq.sigma, rho=propsLiq.rho, w_sup=w_sup_l, eps_liq=eps_liq,
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
     SI.DiffusionCoefficient D_liq[n,nSL]=diffCoeffLiq.D_diluted;

       MediumVapour.DiffusionCoefficient[n] diffCoeffVap(T=T_v, p=p_v[1:n]);
     SI.DiffusionCoefficient D_vap[n,nSV]=diffCoeffVap.D_diluted;
     SI.Velocity w_sup_l[n] "superficial liquid velocity";
     SI.Velocity w_sup_v[n] "superficial vapour velocity";
equation

   k_l = massTransferCoeffLiq.k;
   k_v = massTransferCoeffVap.k;

   alphaLiq=heatTransferCoeffLiq.alpha;
   alphaVap=heatTransferCoeffVap.alpha;

  for j in 1:n loop
     D_av_vap[j]=sum(D_vap[j,:])/nSV;
      D_av_liq[j]=sum(D_liq[j,:])/nSL;
    k_av_vap[j]=sum(massTransferCoeffVap.k[j,:])/nSV;
      k_av_liq[j]=sum(massTransferCoeffLiq.k[j,:])/nSL;
     A_I[j] = interfacialArea.a[j]*geometry.H*geometry.A/n;
  end for;

     /*** superficial velocity ***/
         for j in 1:n loop
  w_sup_l[j] = max(1e-10,Vdot_l[j]/geometry.A);
    w_sup_v[j] = max(1e-10,Vdot_v[j]/geometry.A);
      end for;
end EffectiveDiff;
