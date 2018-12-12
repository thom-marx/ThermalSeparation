within ThermalSeparation.Components.Columns.BaseClasses.Initialization;
model Init_T_xv_p_Ndot0
  "init using T_v, T_l, x_v, x_l_inert, p and Ndot_transfer=0"
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
          T_v = T_v_start;
          T_l = T_l_start;

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
       sum(x_l[j,:])=1;
       sum(x_v[j,:])=1;
    end for;
     T_v = T_v_start;
     T_l = T_l_start;
    p_v[1:n] = p_v_start;
    for i in 1:nSV loop
      if not inertVapour[i] then
        Ndot_v_transfer[:,i]=zeros(n);
      end if;
      end for;

      end if;
end Init_T_xv_p_Ndot0;
