
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<< +++++++++++++++ >>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<< +++++++++++++++ >>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%
%%                                                                            %%
%%                          Dahlquist Test Problem                            %%
%%                                                                            %%
%%  A file that contains all functions necessary for the dahlquist related    %% 
%%  scripts. The simple statement below just prevents Octave from loading the %%
%%  first function definition and immediately executing it.                   %%
1;
%%                                                                            %%
%%                               FUNCTION FILE                                %%
%%                                                                            %%
%%                           Ocatve Script  2.1.0                             %%
%%                                                                            %%
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<< +++++++++++++++ >>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<< +++++++++++++++ >>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%

%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<< IVP Functions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%

    function [rhs] = dahlquist(u)     %% Explicit rhs call
      global lambda;      
      rhs = -lambda*u;
    endfunction

    f = @(x) dahlquist(x);

    function [u] = analytic(t,u0)     %% Provides the analytic solution
      global lambda;    
      u = exp(-lambda*t)*u0;    
    endfunction   
  
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< OSM >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%

    function [un1] = rk1(nt,dt,f,un)
        % Explicit Euler Method -- order 1
        u = un;      
        for i=1:nt
          u = u + dt*feval(f,u);
        endfor
        un1=u;    
    endfunction

    function [un1] = rk3(nt,dt,f,un)
        % Explicit Runge-Kutta Method -- order 3                
        u = un;
        for i=1:nt            
            k1 = feval(f,u);
            k2 = feval(f,u+dt*0.5*k1);
            k3 = feval(f,u-dt*k1+2.0*dt*k2);
            
            u = u +  dt*(k1+4*k2+k3)/6;                  
        endfor
        un1 = u;    
    endfunction

    function [un1] = rk4(nt,dt,f,un)
        % Explicit Runge-Kutta Method -- order 4                
        u = un;
        for i=1:nt
            k1 = feval(f,u);
            k2 = feval(f,u+dt*0.5*k1);
            k3 = feval(f,u+.5*dt*k2);
            k4 = feval(f,u+dt*k3);
            
            u = u +  dt*(k1+2*k2+2*k3+k4)/6;        
        endfor
        un1 = u;    
    endfunction

    function [un1] = cn(nt,dt,un)
      % Crank-Nicolson -- order 2
      global lambda;
      u = un;
      for i=1:nt
        u = (1.0-.5*dt*lambda)*u/(1.0+.5*dt*lambda);
      endfor
      un1 = u;      
    endfunction
    
    function [un1] = irk1(nt,dt,un)
      % Implicit Euler Method -- order 1
        global lambda;
        u = un;
        for i=1:nt
          u = u/(1.0+dt*lambda);
        endfor
        un1 = u;
    endfunction
      
%% <<<<<<<<<<<<<<<<<<<<< Define generic ODE function >>>>>>>>>>>>>>>>>>>>>>>> %%

    function [un1] = call_osm(nt,dt,f,un,method)
        %% Allows for switching comfortably bewtween the OSMs 
        if (method == 0) 
            un1 = rk1(nt,dt,f,un);
        elseif  ( method == 1)
            un1 = rk3(nt,dt,f,un);
        elseif  ( method == 2)
            un1 = rk4(nt,dt,f,un);   
        elseif  ( method == 3)
            un1 = cn(nt,dt,un);
        elseif  ( method == 4)   
            un1 = irk1(nt,dt,un);               
        elseif  ( method == 5)   
            un1 = analytic(nt*dt,un);        
        else
            disp("OSM unkown ... \n");
            abort;
        endif
    endfunction    
    
%% <<<<<<<<<<<<<<<<<<<<<<<<<<< Error Estimates >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%
    
    function [Etot,Ediss,Edisp] = error_estimate(qt,qd)

        M  = max(size(qt));
        v1 = ones(M,1);    
      
        erw_qt   = (1.0/M)*sum(qt);
        erw_qd   = (1.0/M)*sum(qd);
        erw_qtqd = (1.0/M)*sum(qt .* qd );
 
        dev_qt = sqrt( (1.0/M)*sum((qt - erw_qt .* v1) .* (qt - erw_qt .* v1)));
        dev_qd = sqrt( (1.0/M)*sum((qd - erw_qd .* v1) .* (qd - erw_qd .* v1)));      
        
        Etot = (1.0/M)*sum((qt-qd) .* (qt-qd));
        Ediss =  (dev_qt - dev_qd)**2 + (erw_qt - erw_qd)**2;
        Edisp = 2.0*(dev_qt*dev_qd - erw_qtqd + erw_qt*erw_qd);
        
    endfunction
 
%% <<<<<<<<<<<<<<<<<<<<<<< Stability Functions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%

    function [rho] = stabFunc(z,n,method)
      
      zn = z/n;
      if      ( method == 0)
        rhon = 1 + zn;        
      elseif  ( method == 1)
        b    = [1/6 4/6 1/6];
        A    = [0 0 0; .5 0 0; -1 2 0];
                
        rhon = 1 + zn*b*(eye(3,3)-zn*A)*ones(3,1);
      elseif  ( method == 2)
        b = [1/6 1/3 1/3 1/6];
        A = [ 0 0 0 0; 0 .5 0 0; 0 0 .5 0; 0 0 0 1];
        
        rhon = 1 + zn*b*(eye(4,4)+zn*A)*ones(4,1);
      elseif  ( method == 3)
        rhon = ( 1 + zn )/( 1 - zn );
      elseif  ( method == 4)
        rhon = 1/( 1 - zn );
      elseif  ( method == 5)
        rhon = exp(zn);
      else
        disp("OSM unknown ... \n");
        abort;
      endif
    rho = rhon^n;
      
    endfunction
    
%% <<<<<<<<<<<<<<<<<<<<<<< Miscellaneous Functions >>>>>>>>>>>>>>>>>>>>>>>>>> %%

    function [name] = get_osm(method)
      %%
      if      ( method == 0)
        name = "Explicit-Euler";        
      elseif  ( method == 1)
        name = "Runge-Kutta3";
      elseif  ( method == 2)
        name = "Runge-Kutta4";
      elseif  ( method == 3)
        name = "Crank-Nicolson";
      elseif  ( method == 4)
        name = "Implicit-Euler";
      elseif  ( method == 5)
        name = "Exponential-Integrator";        
      else
        disp("OSM unknown ... \n");
        abort;
      endif
    endfunction
    
