within ThermalSeparation.FilmModel.PackedColumn;
model MS "Maxwell-Stefan mass transfer - no film reaction"
extends BaseFilmPacked;
extends ThermalSeparation.FilmModel.BaseClasses.MaxwellStefan;

 final parameter Integer aux[  :] = {1,3,6,10,15, 21, 28, 36, 45};
Geometry geometry(n=n);
replaceable model MassTransferLiquid =
      ThermalSeparation.HeatAndMassTransfer.PackedColumn.Liquid.RochaStructured
                                                                                            constrainedby
    ThermalSeparation.FilmModel.BaseClasses.PackedColumn.BaseMSLiquid        annotation(choicesAllMatching=true,Dialog( tab="Heat and Mass Transfer", group="Liquid and vapour mass transfer"));

      replaceable model MassTransferVapour =
      ThermalSeparation.HeatAndMassTransfer.PackedColumn.Vapour.RochaStructured
        constrainedby
    ThermalSeparation.FilmModel.BaseClasses.PackedColumn.BaseMSVapour
        annotation (
      choicesAllMatching=true, Dialog(
      tab="Heat and Mass Transfer",
      group="Liquid and vapour mass transfer"));

      replaceable model InterfacialArea =
      ThermalSeparation.HeatAndMassTransfer.InterfacialArea.PackedColumn.RochaStructured
                                                                                                        constrainedby
    ThermalSeparation.HeatAndMassTransfer.InterfacialArea.PackedColumn.BasePacked
                                                                                                      annotation(choicesAllMatching=true,Dialog( tab="Heat and Mass Transfer", group="Contact area between liquid and vapour phase"));
    InterfacialArea interfacialArea(
      redeclare record Geometry=Geometry,
      final n=n,rho=propsLiq.rho, sigma=propsLiq.sigma, eta=propsLiq.eta, Vdot_l=Vdot_l, eps_liq=eps_liq);

replaceable model MassTransferCoeffVap =
      ThermalSeparation.HeatAndMassTransfer.PackedColumn.Vapour.RochaStructured
   constrainedby
    ThermalSeparation.HeatAndMassTransfer.PackedColumn.Vapour.BaseVapMT
     annotation(choicesAllMatching=true,Dialog( tab="Heat and Mass Transfer", group="Liquid and vapour mass transfer"));

  ThermalSeparation.FilmModel.BaseClasses.PackedColumn.BaseMSVapour
    massTransferVapour(
     redeclare model MassTransferCoeff=MassTransferCoeffVap,
      redeclare replaceable package MediumVapour =   MediumVapour,
      redeclare record Geometry =   Geometry,
      final n=n, final nS=nSV,  Vdot_v=Vdot_v, Vdot_l=Vdot_l, eta=propsVap.eta, rho=propsVap.rho, x=x_v,
      T= T_v, p=p_v[1:n], sigma=propsLiq.sigma, eps_liq=eps_liq);

replaceable model MassTransferCoeffLiq =
      ThermalSeparation.HeatAndMassTransfer.PackedColumn.Liquid.RochaStructured
   constrainedby
    ThermalSeparation.HeatAndMassTransfer.PackedColumn.Liquid.BaseLiqMT
     annotation(choicesAllMatching=true,Dialog( tab="Heat and Mass Transfer", group="Liquid and vapour mass transfer"));

   ThermalSeparation.FilmModel.BaseClasses.PackedColumn.BaseMSLiquid
    massTransferLiquid(
      redeclare model MassTransferCoeff=MassTransferCoeffLiq,
       redeclare replaceable package MediumLiquid =  MediumLiquid,
       redeclare record Geometry =   Geometry,
       final n=n, final nS=nSL,  Vdot_v=Vdot_v, Vdot_l=Vdot_l, eta=propsLiq.eta, rho=propsLiq.rho,
       x=propsLiq.x, c_star=c_l_star, T = propsLiq.T, p=p_v[1:n], sigma=propsLiq.sigma, eta_comp = eta_comp, eps_liq=eps_liq);

   replaceable model HeatTransferVapour =
      ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient.Vapour.Constant
     constrainedby
    ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient.Vapour.BaseVapour
                                                                                                        annotation(choicesAllMatching=true,Dialog( tab="Heat and Mass Transfer", group="Heat transfer between liquid and vapour phase"));

   HeatTransferVapour heatTransferCoeffVap(
   redeclare replaceable package Medium =  MediumVapour,
   n=n, nS=nSV,  k_av=k_av_vap, D_av=D_av_vap,
   rho=propsVap.rho,  c=c_v, props=propsVap);

     replaceable model HeatTransferLiquid =
      ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient.Liquid.Constant
     constrainedby
    ThermalSeparation.HeatAndMassTransfer.InterfacialHeatTransferCoefficient.Liquid.BaseLiquid
                                                                                                        annotation(choicesAllMatching=true,Dialog( tab="Heat and Mass Transfer", group="Heat transfer between liquid and vapour phase"));

   HeatTransferLiquid heatTransferCoeffLiq(
    redeclare replaceable package Medium = MediumLiquid,
    n=n, nS=nSL, k_av=k_av_liq, D_av=D_av_liq,
    rho=propsLiq.rho,   c=c_l, props=propsLiq);
ThermalSeparation.Units.CoefficentOfMassTransfer k_av_liq[n];
ThermalSeparation.Units.CoefficentOfMassTransfer k_av_vap[n];
SI.DiffusionCoefficient D_av_liq[n];
SI.DiffusionCoefficient D_av_vap[n];
equation
  R_dash_v=massTransferVapour.R_dash;
  R_dash_l=massTransferLiquid.R_dash;
   matrix_beta_v= massTransferVapour.R;
   matrix_beta_l = massTransferLiquid.R;
   alphaLiq=heatTransferCoeffLiq.alpha;
   alphaVap=heatTransferCoeffVap.alpha;

  for j in 1:n loop
     D_av_vap[j]=sum(massTransferVapour.D[j,:])/aux[nSV-1];
      D_av_liq[j]=sum(massTransferLiquid.D[j,:])/aux[nSL-1];
    k_av_vap[j]=sum(massTransferVapour.k_2[j,:])/aux[nSV-1];
      k_av_liq[j]=sum(massTransferLiquid.k_2[j,:])/aux[nSL-1];
     A_I[j] = interfacialArea.a[j]*geometry.H*geometry.A/n;
  end for;
end MS;
