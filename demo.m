%function demo()
    clear all;
    clc;
    
    figure;
    hold;
        
    % initial condition:
    y0 = 0;
    t_start = -0.5;
    t_end = 6;
    p = t_start:0.01:t_end;
 
    
    % Use ode45, reference  
    [tode23t,xode23t]=ode23t(@eqns, p, y0);
    plot(tode23t,xode23t(:,1), '.','color','red','MarkerSize',5);
    
    display ('done ode23t')
    
    % Use Trapsol
    [t,y]=trapsol(@eqns, [t_start t_end],y0);
     plot(t, y,'-','color','black','linewidth',1.4) 
    
    
    legend('ode23t', 'Trapsol');
    xlabel('Time (s)')
   
    
%end