function [J f0] = jacob_fn(jacob_function,ode,t,y,h,previous_values,S1)
            %% Calculates JACOBIAN MATRIX %%%%
            %%% previous value = exact value already calculated for previous time step
            %%% newton_function is the function for which JACOBIAN is
            %%% CALCULATED and is formulated in jacob_function.m
            %%% y = [; ; ; ] current point for calculation
            %%% h = step size current
            %%% t = current value of independent variable  //time
          
            n = size(y,1);
            J = zeros(n,n);
            ee = 1e-6;
            f0 = feval(jacob_function,ode,t,y,h,previous_values,S1);
                                 
            for i = 1:n   
                yj = y;                
                if abs(y(i)) > 1 
                e = ee*abs(y(i));
                else
                    e = ee;
                end
                
                yj(i,1)=y(i,1)+e;
%               S = feval(ode,t,yj);
                S = feval(jacob_function,ode,t,yj,h,previous_values,S1);
                J(:,i) = (S-f0)/e;         
            end
              
end           