within ThermalSeparation.FilmModel.BaseClasses;
package StateSelection

  package StateSelectionNoneq
    "state selection models for nonequilibrium film models"
    partial model BaseStateSelectionNoneq
      "base class for state selection models for nonequilibrium film models"

      final parameter Integer nSL = MediumLiquid.nSubstance;
      final parameter Integer nSV = MediumVapour.nSubstance;
      parameter Integer n(min=1) annotation(Dialog(enable=false));

      replaceable package MediumLiquid = Media.BaseMediumLiquid annotation(Dialog(enable=false));
      replaceable package MediumVapour = Media.BaseMediumVapour annotation(Dialog(enable=false));

      input MediumLiquid.ThermodynamicProperties[n] propsLiq annotation(HideResult=true);
      input MediumVapour.ThermodynamicProperties[n] propsVap annotation(HideResult=true);
      input SI.Pressure p_v[n] annotation(HideResult=true);
      input SI.Concentration c_l[n,nSL] annotation(HideResult=true);

    end BaseStateSelectionNoneq;

    model StateSelection1 "states: c_v, c_l, u_v, u_l, T_v"
      extends ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionNoneq.BaseStateSelectionNoneq;
       SI.Concentration c_v_state[n,nSV](each stateSelect=StateSelect.prefer)= propsVap.c;
       SI.Concentration c_l_state[n,nSL](each stateSelect=StateSelect.prefer) = c_l;
       SI.Concentration u_v_state[n](each stateSelect=StateSelect.prefer) = propsVap.u;
       SI.Concentration u_l_state[n](each stateSelect=StateSelect.prefer) = propsLiq.u;
       SI.Concentration T_v_state[n](each stateSelect=StateSelect.prefer) = propsVap.T;

    end StateSelection1;

    model StateSelection2 "states: c_v, c_l, u_v, u_l, MM_v"
      extends ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionNoneq.BaseStateSelectionNoneq;
      SI.Concentration c_v_state[n,nSV](each stateSelect=StateSelect.prefer) = propsVap.c;
      SI.Concentration c_l_state[n,nSL](each stateSelect=StateSelect.prefer) = c_l;
      SI.Concentration u_v_state[n](each stateSelect=StateSelect.prefer) = propsVap.u;
      SI.Concentration u_l_state[n](each stateSelect=StateSelect.prefer) = propsLiq.u;
      SI.Concentration MM_v_state[n](each stateSelect=StateSelect.prefer) = propsVap.MM;
    equation

    end StateSelection2;

    model None "none"
      extends ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionNoneq.BaseStateSelectionNoneq;

    end None;

    model ReactionEquilibrium "states: c_v, c_l[1:nSL-1], u_v, u_l, MM_v"
      extends ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionNoneq.BaseStateSelectionNoneq;
      SI.Concentration c_v_state[n,nSV](each stateSelect=StateSelect.prefer) = propsVap.c;
      SI.Concentration c_l_state[n,nSL-1](each stateSelect=StateSelect.prefer) = c_l[:,1:nSL-1];
      SI.Concentration u_v_state[n](each stateSelect=StateSelect.prefer) = propsVap.u;
      SI.Concentration u_l_state[n](each stateSelect=StateSelect.prefer) = propsLiq.u;
      SI.Concentration MM_v_state[n](each stateSelect=StateSelect.prefer) = propsVap.MM;
    equation

    end ReactionEquilibrium;

    model StateSelection3 "states: x_v, x_l, p_v, T_l, T_v"
      extends ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionNoneq.BaseStateSelectionNoneq;
       SI.Concentration x_v_state[n,nSV](each stateSelect=StateSelect.prefer) = propsVap.x;
       SI.Concentration x_l_state[n,nSL](each stateSelect=StateSelect.prefer) = propsLiq.x;
       SI.Concentration p_v_state[n](each stateSelect=StateSelect.prefer) = p_v[1:n];
       SI.Concentration T_l_state[n](each stateSelect=StateSelect.prefer) = propsLiq.T;
       SI.Concentration T_v_state[n](each stateSelect=StateSelect.prefer) = propsVap.T;

     //  Real c_v_state[n](each stateSelect=StateSelect.prefer);

    equation
    for j in 1:n loop
     // c_v_state[j] = propsVap[j].c[nSV];
     // x_v_state[j,:] = propsVap[j].x[1:nSV];
      end for;
    end StateSelection3;
  end StateSelectionNoneq;

  package StateSelectionEq
    "state selection models for true equilibrium film models (reduced number of states)"
    partial model BaseStateSelectionEq
      "base class for state selection models for true equilibrium film model (reduced number of states)"

      final parameter Integer nSL = MediumLiquid.nSubstance;
      final parameter Integer nSV = MediumVapour.nSubstance;
      parameter Integer n(min=1) annotation(Dialog(enable=false));

      replaceable package MediumLiquid = Media.BaseMediumLiquid annotation(Dialog(enable=false));
      replaceable package MediumVapour = Media.BaseMediumVapour annotation(Dialog(enable=false));

      input MediumLiquid.ThermodynamicProperties[n] propsLiq;
      input MediumVapour.ThermodynamicProperties[n] propsVap;
      input SI.Pressure p_v[n];
      input SI.Concentration c_l[n,nSL];
    equation

    end BaseStateSelectionEq;

    model SaturatorAbsober "SprayAbsober"
      extends ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionEq.BaseStateSelectionEq;
       SI.Concentration u_v_state[n](each stateSelect=StateSelect.prefer) = propsVap.u;
       SI.Concentration u_l_state[n](each stateSelect=StateSelect.prefer) = propsLiq.u;
      SI.Pressure p_v_state[n](each stateSelect=StateSelect.prefer)= p_v;
     SI.Concentration c_l_state_1[n](each stateSelect=StateSelect.prefer) = c_l[:,2];
    // SI.Concentration c_l_state_2[n](each stateSelect=StateSelect.prefer) = c_l[:,3];
        SI.Concentration T_v_state[n](each stateSelect=StateSelect.prefer) = propsVap.T;

    end SaturatorAbsober;

    model SprayAbsober "SprayAbsober"
      extends ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionEq.BaseStateSelectionEq;
       SI.Concentration u_v_state[n](each stateSelect=StateSelect.prefer) = propsVap.u;
       SI.Concentration u_l_state[n](each stateSelect=StateSelect.prefer) = propsLiq.u;
      SI.Pressure p_v_state[n](each stateSelect=StateSelect.prefer)= p_v;
     SI.Concentration c_l_state_3[n](each stateSelect=StateSelect.prefer) = c_l[:,2];
        SI.Concentration T_v_state[n](each stateSelect=StateSelect.prefer) = propsVap.T;

    end SprayAbsober;

    model StateSelection3 "states: u_v, c_l, p_v "
      extends ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionEq.BaseStateSelectionEq;
      // SI.Concentration c_l_state[n,nSL](stateSelect=StateSelect.prefer) = c_l;
       SI.Concentration u_v_state[n](each stateSelect=StateSelect.prefer) = propsVap.u;
       SI.Concentration u_l_state[n](each stateSelect=StateSelect.prefer) = propsLiq.u;
      SI.Pressure p_v_state[n](each stateSelect=StateSelect.prefer)= p_v;

     //  SI.Concentration T_v_state[n](each stateSelect=StateSelect.prefer) = propsVap.T;

    end StateSelection3;

    model StateSelection4 "states: c_l, T_v , u_l"
      extends ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionEq.BaseStateSelectionEq;
       SI.Concentration c_l_state[n,nSL](each stateSelect=StateSelect.prefer) = c_l;
       SI.Concentration u_l_state[n](each stateSelect=StateSelect.prefer) = propsLiq.u;

          SI.Concentration T_v_state[n](each stateSelect=StateSelect.prefer) = propsVap.T;

    end StateSelection4;

    model None "none"
      extends ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionEq.BaseStateSelectionEq;

    end None;

    model MA "MA"
      extends ThermalSeparation.FilmModel.BaseClasses.StateSelection.StateSelectionEq.BaseStateSelectionEq;
       SI.Concentration u_v_state[n](each stateSelect=StateSelect.prefer) = propsVap.u;
       SI.Concentration u_l_state[n](each stateSelect=StateSelect.prefer) = propsLiq.u;
      SI.Pressure p_v_state[n](each stateSelect=StateSelect.prefer)= p_v;
     SI.Concentration c_l_state_1[n](each stateSelect=StateSelect.prefer) = c_l[:,2];
     SI.Concentration c_l_state_2[n](each stateSelect=StateSelect.prefer) = c_l[:,3];
      SI.Concentration c_l_state_3[n](each stateSelect=StateSelect.prefer) = c_l[:,4];

    end MA;
  end StateSelectionEq;
end StateSelection;
