within ThermalSeparation.Components.Columns.BaseClasses.Initialization;
model Init_x_T_p "init using T_v, T_l, x_v, x_l and p"
extends BaseInit;

initial equation
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
end Init_x_T_p;
