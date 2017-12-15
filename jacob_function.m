function [nfns] =  jacob_function(ode,t,y,h,previous_values,S1)  
%% Formulates f() = 0, i.e. Function for Jacobian Matrix %%%%
% " f() = 0 "  Required form in newton rapshon
%%
%%% Using Euler
%nfns = y-previous_values-h*feval(ode,t,y);    
%%%

%%% Using Trapezoidal
SN = feval(ode,t+h,y);
%S1 = feval(ode,t,previous_values);
nfns = y-previous_values-(h/2)...
    *(SN+S1);
%%%
end
