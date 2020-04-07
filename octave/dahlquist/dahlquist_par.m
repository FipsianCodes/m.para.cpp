
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<< +++++++++++++++ >>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<< +++++++++++++++ >>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%
%%                                                                            %%
%%                          Dahlquist Test Problem                            %%
%%                                                                            %%
%%      Linear Initival Value Problem (IVP) for the Eigenvalue lambda in      %%
%%                             the complex plane                              %%
%%                                                                            %%
%%                            du/dt + lambda*u = 0                            %%
%%                               u(t=0) = u0                                  %%
%%                                                                            %%
%%                         STANDARD PARAREAL ALGORITHM                        %%
%%                                                                            %%
%%                           Ocatve Script  2.1.1                             %%
%%                                                                            %%
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<< +++++++++++++++ >>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<< +++++++++++++++ >>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%

    clear all;

    %% !!! CHANGE TO YOUR WORKING DIRECTORY !!!
    path     = '/home/fips/Dropbox/Dissertation/octave/dahlquist/'; 
    
    src_path = strcat(path,'dahlquist_functions.m');  %% Functions file location
    dat_path = strcat(path,'data/');        %% Data Location
    source(src_path);                       %% Access functions from .m file
    
%% <<<<<<<<<<<<<<<<<<<<<<<<<<< Simulation Setup >>>>>>>>>>>>>>>>>>>>>>>>>>>>> %% 
    
    eps = 1e-12;        %% Stopping criterion Parareal
    
    dT  = .1;           %% Time slice size for Parareal
    nT  = 100;          %% Number of time slices
    NT  = nT+1;         %% Array (vector) size, including the initial condition
    
    ntg = 1;            %% Coarse Solver: Number of sub-timesteps in dT
    nt  = 10;           %% Fine Solver  : Number of sub-timesteps in dT  

    K   = nT;           %% Iteration limit for Parareal
    
    T   = dT*nT;        %% Time Interval T of the IVP
    dtc = dT/ntg;       %% Coarse Solver: sub-timestep size
    dt  = dT/nt;        %% Fine Solver  : sub-timestep size

    %% Output Handler 
    
    save_data = 1;      %% 1 :: save to dat_path -- 0 :: do not safe
    plot_data = 1;      %% 1 :: plot data        -- 0 :: window output only
    
    
%% <<<<<<<<<<<<<<<<<<<<<<< Define Eigenvalue of the IVP >>>>>>>>>>>>>>>>>>>>> %%
    
    global lambda = 1.0*sqrt(-1);
    lr = real(lambda);
    li = imag(lambda);
    
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Initial Value >>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%    

    u0 = 1.0;

%% <<<<<<<<<<<<<<< Define Coarse and Fine Integration Method >>>>>>>>>>>>>>>> %%
    
    mfine = 5;  %%      0 :: Explicit Euler Method
    mcoar = 4;  %%      1 :: Runge-Kutta 3
                %%      2 :: Runge-Kutta 4
                %%      3 :: Crank-Nicolson
                %%      4 :: Implicit Euler Method
                %%      5 :: Exponential Integrator
                
    printf("# ----------------------------------------------------------- #\n");                
    printf("# ----------------------------------------------------------- #\n");                
    printf("#                  Parareal Simulation Setup                  #\n"); 
    printf("# ----------------------------------------------------------- #\n");                    
    printf("# ----------------------------------------------------------- #\n");                    
    printf("#                                                              \n");    
    printf("# Time Interval T of the IVP  :  %2.3f                       \n",T);                
    printf("# Intival Value u0            :  %4.2f                      \n",u0);                
    printf("# Number of time slices       :  %d                         \n",nT);                
    printf("# Time slice size             :  %2.3f                      \n",dT);  
    printf("# Coarse: sub-timestep size   :  %2.3f                     \n",dtc);  
    printf("# Fine  : sub-timestep size   :  %2.3f                      \n",dt);  
    printf("# Parareal stopping criterion :  %1.1e                     \n",eps);  
    printf("# Eigenvalue lambda  (real)   :  %1.2f + %1.2fi          \n",lr,li);    
    printf("# Coarse Solver               :  %s             \n",get_osm(mcoar));    
    printf("# Fine   Solver               :  %s             \n",get_osm(mfine));    
    printf("#                                                              \n");       
    printf("# ----------------------------------------------------------- #\n");      
  
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< +++++++++++ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Run Program >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< +++++++++++ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%
    
%% <<<<<<<<<<<<<<<<<<<<<<<< Reference Solution >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%

    uref = zeros(NT,2);
    uref(1,2) = u0;
    
    for i=1:nT
      uref(i+1,1) = dT*i;
      uref(i+1,2) = call_osm(nt,dT/nt,f,uref(i,2),mfine);
    endfor
    
    if (save_data == 1)
      urefi = [uref(:,1),imag(uref(:,2))];
      urefr = [uref(:,1),real(uref(:,2))];
      save([path,'Uref_i.dat'],'urefi','-ascii'); %% Real part 
      save([path,'Uref_r.dat'],'urefr','-ascii'); %% Imag part
    endif     

%% <<<<<<<<<<<<<<<<<<<<<<<<<<< Fine Solution >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%

    ufine = zeros(NT,2);
    ufine(1,2) = u0;

    for i=1:nT
      ufine(i+1,1) = dT*i;
      ufine(i+1,2) = call_osm(1,dT,f,ufine(i,2),5);
    endfor
 
    if (save_data == 1) 
      ufinei = [ufine(:,1),imag(ufine(:,2))];       
      ufiner = [ufine(:,1),real(ufine(:,2))];
      save([path,'Ufine_i.dat'],'ufinei','-ascii'); %% Real part 
      save([path,'Ufine_r.dat'],'ufiner','-ascii'); %% Imag part
    endif
           
    printf("# --------------------- Start Parareal ---------------------- #\n");          
    printf("# ----------------------------------------------------------- #\n"); 
    
%% <<<<<<<<<<<<<<<<<<<<<<<<<<< First Iteration >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%

    err    = zeros(K,4);    %% Allocate array for convergence evaluation
    Uiter  = zeros(NT,K);   %%  
    Uhat   = zeros(NT,K);   %% Fine Solver, parallel computation
    Utilde = zeros(NT,K);   %% Coarse Solver, serial computation 
    
    for k=1:K
      Uiter(1,k)  = u0;     %% Assign u0 for all iterations K
      Uhat(1,k)   = u0;
      Utilde(1,k) = u0;    
    endfor
 
    for n = 1:nT            %% Assign initial condition for all time slices 
      Utilde(n+1,1) = call_osm(ntg,dT/ntg,f,Utilde(n,1),mcoar); 
      Uiter(n+1,1) = Utilde(n+1,1);
    endfor
    
    err(1,1) = 0;
    err(1,2) = abs(norm(ufine(:,2)- Uiter(:,1),"inf")); %% Error to analytic
    err(1,3) = abs(norm(uref(:,2) - Uiter(:,1),"inf")); %% Error to Reference
    err(1,4) = abs(norm(Uiter(:,1),"inf"));             %% Iteration error
    printf('# It.No. = %i  dlt = %1.1e   e_a = %1.1e   e_r = %1.1e  \n',...
           0,err(1,4),err(1,2),err(1,3));
    
%% <<<<<<<<<<<<<<<<<<<<<<<<<< PARAREAL >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%
           
    for k = 1:K-1

    %% Fine Computation (in Parallel) of F(k)
    
        for n = k:nT
          Uhat(n+1,k) = call_osm(nt,dt,f,Uiter(n,k),mfine);
        endfor
    
        Uiter(k+1,k:K) = Uhat(k+1,k);
    
    %% Iteration procedure          
    
        for n = k:nT
          %% Serial computation of G(k+1) 
          Utilde(n+1,k+1) = call_osm(ntg,dT/ntg,f,Uiter(n,k+1),mcoar);
          %% Update U = G(k+1) + F(k) - G(k)
          Uiter(n+1,k+1) = Utilde(n+1,k+1) + Uhat(n+1,k) - Utilde(n+1,k);
        end
        
     %% Error estimates 
     
        err(k+1,1) = k;
        err(k+1,2) = abs(norm(ufine(:,2)- Uiter(:,k+1),"inf"));
        err(k+1,3) = abs(norm(uref(:,2) - Uiter(:,k+1),"inf"));
        err(k+1,4) = abs(norm(  Uiter(NT,k+1) - Uiter(NT,k), "inf"));
      
        printf('# It.No. = %i  dlt = %1.1e   e_a = %1.1e   e_r = %1.1e  \n',...
               k,err(k+1,4),err(k+1,2),err(k+1,3));
     
    %% Stopping Criterion
    
        if ( err(k+1,4) < eps  )  
            K_break = k;
            break;
        end
        K_break = k;
    end

    err = err(1:K_break+1,:);

    printf("# ----------------------------------------------------------- #\n");            
    printf("# --------------------- Task completed ---------------------- #\n");          
    printf("# ----------------------------------------------------------- #\n");    
    
    
    if ( save_data == 1)
      lastit = [uref(:,1) , real(Uiter(:,K_break)) , imag(Uiter(:,K_break)) ];
      save([dat_path,'last_iter.dat'],'lastit','-ascii');
      
      errname = strcat('err_',get_osm(mcoar),'_',num2str(real(lambda)),'_',num2str(imag(lambda)));
      save([dat_path,errname],'err','-ascii');
    endif
    
    if ( plot_data == 1)
      figure(1);
      semilogy(err(:,1),err(:,2),'-s',err(:,1),err(:,3),'-o',...
               err(:,1),err(:,4),'-x');
      title("Parareal Convergence Behaviour");   
      xlabel("Iterations k");
      ylabel("Inf-Error");
      legend(["Error to Reference";"Error to Analytic";"Error of Iterations"]);
      
      
      if ( K_break > 4)
        plot_count = idivide(K_break, 4, "fix");
      else
        plot_count = 1;
      endif
      
      figure(2);
      hold on;
      plot(uref(:,1),uref(:,2),'-o');
      for i=1:3
        plot(uref(:,1),Uiter(:,i*plot_count));
      endfor
      plot(uref(:,1),Uiter(:,K_break),'-s');
      hold off;
      title("Solution and Iterations");
      xlabel("Time t");
      ylabel("Real(u)");
      legend(['Uref';strcat('Iter ',num2str(plot_count)); ...
                     strcat('Iter ',num2str(2*plot_count)); ...
                     strcat('Iter ',num2str(3*plot_count));'Last Iter']);
    endif
    
clear;


