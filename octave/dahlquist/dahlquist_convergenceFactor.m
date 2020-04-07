
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
%%                   EVALUATION OF THE CONVERGENCE FACTOR                     %%
%%                                                                            %%
%%                           Ocatve Script  2.1.2                             %%
%%                                                                            %%
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<< +++++++++++++++ >>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<< +++++++++++++++ >>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%

    clear all;

    %% !!! CHANGE TO YOUR WORKING DIRECTORY !!! %%
    path     = '/home/fips/Dropbox/Dissertation/octave/dahlquist/'; 
    
    src_path = strcat(path,'dahlquist_functions.m');  %% Functions file location
    source(src_path);                       %% Access functions from .m file
    
    dat_path = strcat(path,'data_stab/');   %% Data Location
    img_path = strcat(path,'img/');         %% Image Location
    
%% <<<<<<<<<<<<<<<<<<<<<<<<<<< Simulation Setup >>>>>>>>>>>>>>>>>>>>>>>>>>>>> %% 
    
    n  = 151;         %% Field length of real(z)
    m  = 251;         %% Field width of imag(z)
    
    nc = 1;           %% Amount of timesteps of the coarse solver 
    nf = 10;          %% Amount of timesteps of the fine solver

    
    lim_real = 5;     %% real(z) in [-lim_real         0]
    lim_imag = 15;    %% imag(z) in [-lim_imag  lim_imag]

%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<< Output Handler >>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%
    
    save_data = 1;      %% 1 :: save to dat_path -- 0 :: do not safe
    plot_data = 1;      %% 1 :: plot data        -- 0 :: window output only 
 
    level_max = 1.1;
    level     = level_max/11;
    contour_levels = [0:level:level_max]; %% Set the levels for the contour plot  
    
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<< Allocate Fields >>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%    

    RHO = zeros(n,m); %% rho(z)
    ZR  = zeros(n,1);   %% real(z)
    ZI  = zeros(m,1);   %% imag(z)
    
%% <<<<<<<<<<<<<<< Define Coarse and Fine Integration Method >>>>>>>>>>>>>>>> %%

    mfine = 5;  %%      0 :: Explicit Euler Method
    mcoar = 0;  %%      1 :: Runge-Kutta 3
                %%      2 :: Runge-Kutta 4
                %%      3 :: Crank-Nicolson
                %%      4 :: Implicit Euler Method
                %%      5 :: Exponential Integrator

    printf("# ----------------------------------------------------------- #\n");                
    printf("# ----------------------------------------------------------- #\n");                
    printf("#             Convergence Analysis Simulation Setup           #\n"); 
    printf("# ----------------------------------------------------------- #\n");                    
    printf("# ----------------------------------------------------------- #\n");                    
    printf("#                                                              \n");                
    printf("# Intival Value u0            :  %4.2f                       \n",1);  
    printf("# Coarse: time steps          :  %d                         \n",nc);  
    printf("# Fine  : time steps          :  %d                         \n",nf);    
    printf("# Coarse Solver               :  %s             \n",get_osm(mcoar));    
    printf("# Fine   Solver               :  %s             \n",get_osm(mfine));    
    printf("#                                                              \n");       
    printf("# ----------------------------------------------------------- #\n");                 
                
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< +++++++++++ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Run Program >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%
%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< +++++++++++ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %%                
    
    %%% Generating the z-field for the upcoming computation of the rho(z) field
        
    for i=1:n
        lambda_real = -lim_real + (i-1)*lim_real/(n);
        ZR(i) = lambda_real;      
    endfor
    for j=1:m
        lambda_cplx = -lim_imag + 2.*(j-1)*lim_imag/(m-1);
        ZI(j) = lambda_cplx;
    endfor
    ZZR = ZR * (ones(1,m));
    ZZI = ones(n,1) * ZI';
    Z = ZZR + sqrt(-1)*ZZI;     %% z-field of size nxm
    
  
    %%% Computing the convergence factor rho(z)
                
    for i=1:n
      
      for j=1:m
          
          Rf = stabFunc(Z(i,j),nf,mfine);

          Rc = stabFunc(Z(i,j),nc,mcoar);
          
          rho = abs(Rf - Rc)/(1 - abs(Rc));
          RHO(i,j) = rho;
             
          ZI(j) = lambda_cplx;            
      endfor
               
    endfor            
    
    %%% Data Plotting 
    if ( plot_data == 1 )
      conv_plot = figure(1);    
      colormap ("jet");
      [c, h] = contour(RHO',contour_levels);
      clabel (c, h, contour_levels, "fontsize", 12);
      set (gca, 'xtick', [1 n/2 n]);
      set (gca, 'xticklabel', {'-5', '-2.5', '0'});
      set (gca, 'ytick', [1 round(m/4) round(m/2) round(m*3/4) m]);
      set (gca, 'yticklabel', {'-15', '-7.5', '0', '7.5', '15'});
      xlabel('Re(z)');
      ylabel('Im(z)');
      if ( save_data == 1 )
          file_name = strcat("conv_",get_osm(mfine),"_",get_osm(mcoar),".eps");
          print (conv_plot, strcat(img_path,file_name),"-color");
      endif
    endif
    
    printf("# ---------------------- Task complete! --------------------- #\n");
    printf("# ----------------------------------------------------------- #\n");
    
clear;    

