within ThermalSeparation.Components.Columns.BaseClasses;
partial model FeedColumn
  "column with optional liquid and or vapour feed inlets"
 //extends BaseColumn;
 extends ThermalSeparation.Components.Columns.BaseClasses.BaseColumn;
// extends ThermalSeparation.Components.Columns.BaseClasses.BaseColumn_external;
/***LIQUID FEED ***/

//protected
model InternalFeedPort "internal model for conditional stream connectors"
  ThermalSeparation.Interfaces.LiquidPortOut          feedLiquidInternal;
  ThermalSeparation.Interfaces.GasPortOut          feedVapourInternal;
end InternalFeedPort;
  InternalFeedPort internalFeedPort[n](feedLiquidInternal(redeclare each
        package Medium =
          MediumLiquid),feedVapourInternal(redeclare each package Medium =
          MediumVapour));
protected
  parameter Integer[n-numberLiquidFeedsInternal] aux1(each fixed=false);
  parameter Integer[n] aux2(each fixed=false);
  parameter Integer counterL(fixed = false)
    "counting the stages; equals n at end of algorithm";

public
  parameter Boolean hasLiquidFeed = false "true, if there exist a liquid feed" annotation(Dialog(tab="Feed", group="Liquid Feed"));
  parameter Integer numberLiquidFeeds(min=0,max=n) = 1  annotation(Dialog(enable = hasLiquidFeed,tab="Feed", group="Liquid Feed"));
  parameter Integer[numberLiquidFeeds] stageLiquidFeed={2}
    "number of stage where feed enters the column" annotation(Dialog(enable = hasLiquidFeed, tab="Feed", group="Liquid Feed"));
  final parameter Integer numberLiquidFeedsInternal(min=0,max=n) = if hasLiquidFeed then numberLiquidFeeds else 0;
  ThermalSeparation.Interfaces.LiquidPortIn[numberLiquidFeeds]
                                                      feedLiquid(redeclare
      each package Medium =
        MediumLiquid) if hasLiquidFeed
                     annotation (Placement(transformation(
          extent={{-94,-14},{-74,6}}, rotation=0), iconTransformation(extent={{
            -94,-14},{-74,6}}))); // ;

  ThermalSeparation.Utilities.MediumLink mediumLink[n];
  MediumLiquid.BaseProperties mediumLiquidFeed[numberLiquidFeeds](each T0=T_ref, each p=p_v[n+1],h=actualStream(feedLiquid_dummy.h_outflow), x=actualStream(feedLiquid_dummy.x_outflow)) if                  hasLiquidFeed;

/***VAPOUR FEED ***/
//   ThermalSeparation.Interfaces.GasPortOut[n] feedVapourInternal(redeclare each
//       package Medium =
//         MediumVapour,                                                                                       x(start=x_v_start), c(start=c_v));
protected
  parameter Integer[n-numberVapourFeedsInternal] testV(each fixed=false);
  parameter Integer[n] test2V(each fixed=false);
  parameter Integer counterLV(fixed = false)
    "counting the stages; equals n at end of algorithm";

public
  parameter Boolean hasVapourFeed = false "true, if there exist a liquid feed" annotation(Dialog(tab="Feed", group="Vapour Feed"));
  parameter Integer numberVapourFeeds(min=0,max=n) = 1  annotation(Dialog(enable = hasVapourFeed, tab="Feed", group="Vapour Feed"));
  parameter Integer[numberVapourFeeds] stageVapourFeed={1}
    "number of stage where feed enters the column" annotation(Dialog(enable = hasVapourFeed,tab="Feed", group="Vapour Feed"));
  final parameter Integer numberVapourFeedsInternal(min=0,max=n) = if hasVapourFeed then numberVapourFeeds else 0;
//   ThermalSeparation.Utilities.LinkVapourSink[n] linkVapour(redeclare each
//       package Medium =
//         MediumVapour,   p=p_v[1:n]);
  ThermalSeparation.Interfaces.GasPortIn[numberVapourFeeds]
                                                   feedVapour(redeclare each
      package Medium =
        MediumVapour) if hasVapourFeed annotation (Placement(transformation(
          extent={{-94,8},{-74,28}}, rotation=0), iconTransformation(extent={{
            -94,8},{-74,28}})));

  ThermalSeparation.Utilities.MediumLink mediumVapourLink[n];
  MediumVapour.BaseProperties mediumVapourFeed[numberVapourFeeds](each T0=T_ref, each p=p_v[n+1],    c=c_v_feed_used,                         h=actualStream(feedVapour_dummy.h_outflow), x=actualStream(feedVapour_dummy.x_outflow),  x_star=actualStream(feedVapour_dummy.x_outflow)) if                   hasVapourFeed;
                                                                                                /*, TODO FIXME! c=c_v_feed */
 // Modelica.Blocks.Interfaces.RealInput h_help;

  ThermalSeparation.Interfaces.GasPortIn[numberVapourFeeds] feedVapour_dummy(redeclare
      each package                                                                                  Medium =
        MediumVapour);
  ThermalSeparation.Interfaces.LiquidPortIn[numberLiquidFeeds] feedLiquid_dummy(redeclare
      each package                                                                                     Medium =
        MediumLiquid);
  ThermalSeparation.Interfaces.GasPortOut[numberVapourFeeds] feedVapour_dummy2(redeclare
      each package                                                                                    Medium =
        MediumVapour) if                                                                                                        not hasVapourFeed;
  ThermalSeparation.Interfaces.LiquidPortOut[numberLiquidFeeds] feedLiquid_dummy2(redeclare
      each package                                                                                       Medium =
        MediumLiquid) if                                                                                                            not hasLiquidFeed;
  SourcesSinks.SourceGas sourceGas[numberVapourFeeds](each use_Flow=false, redeclare
      each package                                                                                Medium =
        MediumVapour,                                                                                                 each T=293.15,each x = x_v_start_const,each Flow=1) if not hasVapourFeed;
  SourcesSinks.SourceLiquid sourceLiquid[numberLiquidFeeds](each use_Flow=false,redeclare
      each package                                                                                     MediumLiquid =
        MediumLiquid,                                                                                                            each T=293.15,each x=x_l_start_const,each Flow=1) if not hasLiquidFeed;
  SourcesSinks.SinkGas sinkGas[numberVapourFeeds](redeclare each package Medium =
        MediumVapour,                                                                        each p=100000) if not hasVapourFeed;
  SourcesSinks.SinkLiquid sinkLiquid[numberLiquidFeeds](redeclare each package Medium =
        MediumLiquid,                                                                               each p=100000) if not hasLiquidFeed;

  Real c_v_feed_used[numberVapourFeeds,nSV];

equation
  for i in 1:numberVapourFeeds loop
    if hasVapourFeed then
      connect(feedVapour_dummy[i], feedVapour[i]);
    else
      connect(feedVapour_dummy[i], sourceGas[i].gasPortOut);
      connect(feedVapour_dummy[i], feedVapour_dummy2[i]);
      connect(feedVapour_dummy2[i],sinkGas[i].gasPortIn);
    end if;
  end for;
  for i in 1:numberLiquidFeeds loop
    if hasLiquidFeed then
      connect(feedLiquid_dummy[i], feedLiquid[i]);
    else
      connect(feedLiquid_dummy[i],sourceLiquid[i].liquidPortOut);
      connect(feedLiquid_dummy[i], feedLiquid_dummy2[i]);
      connect(feedLiquid_dummy2[i],sinkLiquid[i].liquidPortIn);
    end if;
  end for;

/*** link to base class ***/

for j in 1:n loop
Vdot_v_feed[j] = internalFeedPort[j].feedVapourInternal.Ndot/sum(c_v_feed[j,:]);
c_v_feed[j,:] = inStream(internalFeedPort[j].feedVapourInternal.x_outflow[:])*rho_v_feed[j]/MM_v_feed[j];
h_v_feed[j] = inStream(internalFeedPort[j].feedVapourInternal.h_outflow);
rho_v_feed[j] = mediumVapourLink[j].mediumConIn.rho;
MM_v_feed[j] = mediumVapourLink[j].mediumConIn.MM;
Vdot_l_feed[j] = internalFeedPort[j].feedLiquidInternal.Ndot/sum(c_l_feed[j,:]);
c_l_feed[j,:] = inStream(internalFeedPort[j].feedLiquidInternal.x_outflow[:])*rho_l_feed[j]/MM_l_feed[j];
h_l_feed[j] = inStream(internalFeedPort[j].feedLiquidInternal.h_outflow);
rho_l_feed[j] = mediumLink[j].mediumConIn.rho;
MM_l_feed[j] = mediumLink[j].mediumConIn.MM;
end for;

/***LIQUID FEED ***/
//conditional connectors feed
for j in 1:numberLiquidFeeds loop
  connect(feedLiquid[j],internalFeedPort[stageLiquidFeed[j]].feedLiquidInternal);
  connect(mediumLink[stageLiquidFeed[j]].mediumConIn,mediumLiquidFeed[j].mediumConOut);
end for;
for j in 1:n loop
  internalFeedPort[j].feedLiquidInternal.p = p_v[j];
end for;
if hasLiquidFeed then
  for j in 1:(n-numberLiquidFeeds) loop
    internalFeedPort[aux1[j]].feedLiquidInternal.h_outflow = h_l[aux1[j]];//1e5;
    internalFeedPort[aux1[j]].feedLiquidInternal.x_outflow = x_l[aux1[j],:];//1/nSL * ones(nSL);//

    mediumLink[aux1[j]].mediumConIn.h = 1;
    mediumLink[aux1[j]].mediumConIn.rho = 1;
    mediumLink[aux1[j]].mediumConIn.MM = 1;
  end for;
  for j in 1:numberLiquidFeeds loop
    internalFeedPort[stageLiquidFeed[j]].feedLiquidInternal.h_outflow = h_l[stageLiquidFeed[j]];//1e5;
    internalFeedPort[stageLiquidFeed[j]].feedLiquidInternal.x_outflow = x_l[stageLiquidFeed[j],:];//1/nSL * ones(nSL);
  end for;
else
  for j in 1:n loop
    internalFeedPort[j].feedLiquidInternal.h_outflow = h_l[j];//1e5;
    internalFeedPort[j].feedLiquidInternal.x_outflow = x_l[j,:];//nSL * ones(nSL);

    mediumLink[j].mediumConIn.h = 1;
    mediumLink[j].mediumConIn.rho = 1;
    mediumLink[j].mediumConIn.MM = 1;
  end for;
end if;

/***VAPOUR FEED ***/
//conditional connectors feed
for j in 1:numberVapourFeeds loop
  connect(feedVapour[j],internalFeedPort[stageVapourFeed[j]].feedVapourInternal);
  connect(mediumVapourLink[stageVapourFeed[j]].mediumConIn,mediumVapourFeed[j].mediumConOut);
end for;
for j in 1:n loop
  internalFeedPort[j].feedVapourInternal.p = p_v[j];
end for;
if hasVapourFeed then
for j in 1:(n-numberVapourFeeds) loop
  internalFeedPort[testV[j]].feedVapourInternal.h_outflow = h_v[testV[j]];//1e5;//
  internalFeedPort[testV[j]].feedVapourInternal.x_outflow = x_v[testV[j],:];//1/nSV * ones(nSV);//

  mediumVapourLink[testV[j]].mediumConIn.h = 1;
  mediumVapourLink[testV[j]].mediumConIn.rho = 1;
  mediumVapourLink[testV[j]].mediumConIn.MM = 1;
end for;
  for j in 1:numberVapourFeeds loop
    internalFeedPort[stageVapourFeed[j]].feedVapourInternal.h_outflow = h_v[stageVapourFeed[j]];//1e5;
    internalFeedPort[stageVapourFeed[j]].feedVapourInternal.x_outflow = x_v[stageVapourFeed[j],:];//
  end for;
else
  for j in 1:n loop
    internalFeedPort[j].feedVapourInternal.h_outflow = h_v[j];//1e5;
    internalFeedPort[j].feedVapourInternal.x_outflow = x_v[j,:];//nSL * ones(nSL);

    mediumVapourLink[j].mediumConIn.h = 1;
    mediumVapourLink[j].mediumConIn.rho = 1;
    mediumVapourLink[j].mediumConIn.MM = 1;
  end for;
end if;

for i in 1:numberVapourFeeds loop
    c_v_feed_used[i]= c_v_feed[stageVapourFeed[i],:];
end for;

/***LIQUID FEED ***/
initial algorithm
aux2 :={i for i in 1:n};
for j in 1:numberLiquidFeedsInternal loop
    for i in j:n loop
        if stageLiquidFeed[j]==aux2[i] then
            aux2[i]:=0*aux2[i];
        end if;
    end for;
end for;
counterL:=1;
for g in 1:n loop
    if aux2[g]>0 then
        aux1[counterL]:=aux2[g];
        counterL:=counterL + 1;
    end if;
end for;

/***VAPOUR FEED ***/
test2V :={i for i in 1:n};
for j in 1:numberVapourFeedsInternal loop
    for i in j:n loop
        if stageVapourFeed[j]==test2V[i] then
            test2V[i]:=0*test2V[i];
        end if;
    end for;
end for;
counterLV:=1;
for g in 1:n loop
    if test2V[g]>0 then
        testV[counterLV]:=test2V[g];
        counterLV:=counterLV + 1;
    end if;
end for;

  annotation (         Icon(coordinateSystem(preserveAspectRatio=false, extent=
            {{-100,-100},{100,100}}),
                            graphics),
    Documentation(info="<html>
<p>This class provides the equations necessary to describe vapour and / or liquid feeds to the column.</p>
</html>"));
end FeedColumn;
