within ThermalSeparation.Components.Columns.BaseClasses;
package Initialization 
  model Dyncap1_GG
  extends ThermalSeparation.Components.Columns.BaseClasses.Initialization.BaseInit;
  initial equation
     for j in 1:n loop
           for i in 1:nSL-2 loop
                 //x_l[j,i] = x_l_start[j,i];
           end for;
          for i in 2:nSV loop
              x_v[j,i] = x_v_start[j,i];
          end for;
         //sum(x_l[j,:])=1;
         x_l[j,3]/(x_l[j,1]*0.018)=7;
         sum(x_v[j,:])=1;
      end for;
      T_v[1:n] = T_v_start;
      // T_l[1:n] = T_l_start;
      p_v[1:n] = p_v_start;
         for i in 1:nSV loop
           if not inertVapour[i] then
             //Ndot_v_transfer[:,i]=zeros(n);
           end if;
         end for;
  end Dyncap1_GG;

  model DyncapStartUpAbsorption "trueequilibriumStartUpDyncapAbsorber"
    import ThermalSeparation;
  extends BaseInit;
  initial equation
     for j in 1:n loop
           for i in 1:nSL-2 loop
                 x_l[j,i] = x_l_start[j,i];
           end for;
          for i in 1:nSV-1 loop
              x_v[j,i] = x_v_start[j,i];
          end for;
         //sum(x_l[j,:])=1;
         x_l[j,3]/(x_l[j,1]*0.018)=7;
         sum(x_v[j,:])=1;
      end for;
       //T_v[1:n] = T_v_start;
      T_l[1:n] = T_l_start;
      p_v[1:n] = p_v_start;
         for i in 1:nSV loop
           if not inertVapour[i] then
             //Ndot_v_transfer[:,i]=zeros(n);
           end if;
         end for;
  end DyncapStartUpAbsorption;

  model DyncapStartUpDesorption "trueequilibriumStartUpDyncapDesorber"
  extends BaseInit;
  initial equation
          for j in 1:n loop
              for i in 1:nSL-2 loop
              c_l[j,i]=  x_l_start[j,i] /v[j];
               //x_l[j,i] = x_l_start[j,i];
              end for;
              for i in 1:nSV loop
               x_v[j,i] = x_v_start[j,i];
              end for;
                  x_l[j,3]/(x_l[j,1]*0.018)=7;
         sum(x_l[j,:])=1;
         //sum(x_v[j,:])=1;
          end for;
            T_v = T_v_start;
            //T_l = T_l_start;
          // p_hyd[1:n] = p_v_start;
           //oder
           p_v[1:n] = p_v_start;
  end DyncapStartUpDesorption;
end Initialization;
