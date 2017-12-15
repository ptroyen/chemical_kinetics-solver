function [t,y] = trapsol(ode,tspan,y0)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS IS THE ODE SOLVER
% functions for differential equations: ode
% tspan
% initial conditions / y0 = [ , , , ];
%%% Output
% Time
% Solution of y with time: in this form "y(:(time),k(no_of_outputs))"
% Solution of t in this form : "t = [ ; ; ; ; ]"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

previous_values = y0';                      %%% Initial Conditions
t_start = tspan(1,1);
t_end = tspan(1,2);
t = [];
y = [];

% Save Starting values
y(:,1) = previous_values;
t(1,1) = t_start;
kk = 0;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
    h = 1e-2;          %%% STEP SIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%

    th = t_start;
    
    k = 1;
    
    
    while(th<t_end-h)

        %% Initial Guess using RK2
        S1 = feval(ode,th,previous_values);
        S2 = feval(ode,th+h,previous_values+h*S1); 
        new_value = previous_values+ h*(S1+S2)/2;
  
        %% Newton-Raphson-Jacobian Iteration Starts
        while(1)  
            %%% previous value = exact value already calculated for previous time step
            %%% new value = guess value for this time step
            kk = kk+1;
            
            %%% Calculate the Jacobian matrix  
            %SN = feval(ode,t+h,y);
            [J f0] = jacob_fn(@jacob_function,ode,th,new_value,h,previous_values,S1); 
 
            %%% solving for delta 
            delta = - J\f0;
            
            %%% New Value
            check_value = delta + new_value;
            check = abs(delta)/check_value;
             
        
           %%% Terminating Condition  
           if max(check) <1e-7        %% 8 10 12 14 15
                break;
           end
           new_value = check_value;   
           end
    
    
         th = th + h;
         k = k+1;
         SN = feval(ode,th,check_value);
         t(1,k) = th;
         y(:,k) = previous_values+(h/2)*(SN+S1);                                                
         previous_values = y(:,k);
         
         
    end
     y = y';
     t = t';
     
end