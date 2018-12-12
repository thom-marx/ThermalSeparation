within ThermalSeparation.Components.Columns.BaseClasses.Initialization;
model Init_EquilibriumBalancing "for equilibrium balancing model"
extends BaseInit;
initial equation
  if considerStartUp then
        for j in 1:n loop
            for i in 1:nSL loop
           c_l[j,i]=  x_l_start[j,i] /v[j];
            end for;
            for i in 1:nSV loop
             x_v[j,i] = x_v_start[j,i];
            end for;
          end for;
         // T_v = T_v_start;
          T_l = T_l_start;
        // p_hyd[1:n] = p_v_start;
         //oder
         p_v[1:n] = p_v_start;

  else

            for j in 1:n loop
          for i in 1:nSL loop
             if inertLiquid[i] then
                x_l[j,i] = x_l_start[j,i];
             end if;
          end for;
               for i in 3:nSV loop
              x_v[j,i] = x_v_start[j,i];
          end for;
//        sum(x_l[j,:])=1;
//        sum(x_v[j,:])=1;
    end for;
   /*** only one temperature can be initialized ***/
    T_l = T_l_start;
    p_v[1:n] = p_v_start;

      end if;
end Init_EquilibriumBalancing;
