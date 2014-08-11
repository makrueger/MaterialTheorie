%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Axial stress driver for tension-compression-simulations
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% 25.06.2010; tw
% 15.06.2011: rh
% 09.07.2012: rb
% 03.07.2014: rh
%-------------------------------------------------------------------------

clc
clear all
close all

%% set flags to control calculation ---------------------------------------

% 1: Hooke's law
% 2: Constitutive law of Kauderer
% 3: Maxwell model
% 4: Smallstrain von Mises plasticity
matflag = 4;

% calculation of moduli
% 1: analytically
% 2: numerically
mflag = 1;


%% prescribed material parameters -----------------------------------------
if matflag == 1
    emod   = 2.1e5;
    nu     = 0.25;
    mat(1) = (emod*nu)/((1+nu)*(1-2*nu)); % lambda
    mat(2) = emod/(2*(1+nu));             % mu
elseif matflag == 2
    mat(1) = 8.3333e4; % K
    mat(2) = 3.8461e4; % G
    mat(3) = 0.5;      % kapvol
    mat(4) = 0.25;     % kapdev
elseif matflag == 3
    mat(1) = 170.83; % K
    mat(2) = 78.85;  % G
    mat(3) = 100;    % kappa
    mat(4) = 40;     % mu
elseif matflag == 4
    mat(1) = 205000; % Young's modulus
    mat(2) = 0.29;   % Poisson ratio
    mat(3) = 2091;   % Hardening modulus
    mat(4) = 695;    % Initial yield stress
    mat(5) = 0.5;    % r=1: isotropic hardening, r=0: kinematic hardening
else
    error('Choose matflag = 1 ... 4')
end;


%% prescribed loading path ------------------------------------------------
% 1: linear increasing
% 2: linear increasing with plateau afterwards
% 3: cyclic loading (2 cycles)
% 4: cyclic loading (4 cycles)
% 5: positive cyclic loading (4 cycles)

ltype=4;

emax = 0.010; % maximum strain
if ltype==1
    t   =[0    1    7];
    lam =[0    emax  emax];
elseif ltype==2
    t   =[0     1     3    4      6      7     9  15];
    lam =[0  emax  emax 2*emax 2*emax 3*emax 3*emax 3*emax];
elseif ltype==3
    t  = [0   0.5   1.5   2.5   3.5  4];
    lam= [0  emax -emax  emax -emax  0];
elseif ltype==4
    t  = [0   0.5   1.5   2.5   3.5   4.5   5.5   6.5   7.5  8];
    lam= [0  emax -emax  emax -emax  emax -emax  emax -emax  0];
elseif ltype==5
    t  = [0   0.5   1.5   2.5   3.5   4.5   5.5   6.5   7.5  8];
    lam= [0  emax 0  emax 0  emax 0  emax 0  0];
else
    error('Choose ltype = 1 ... 5')
end

% prescribed load/time step
dt = 0.01;

% start and end-time of loading, time-scale, no. of steps
ta=t(1);
te=t(end);
time=ta:dt:te;
steps=size(time,2);
e11=load_steplin(dt,t,lam);


%% set initial values of internal variables -------------------------------
sdv.ep = zeros(3,3);
sdv.a  = zeros(3,3);
sdv.k  = 0;

% initialise strain and stress tensors
e = zeros(3,3);
s = zeros(3,3);

% initialise quantities for post-processing
s11      =zeros(steps,1);
ep11     =zeros(steps,1);
sviscoin = zeros(3,3);

%% main loop --------------------------------------------------------------

% Newton-tolerance and maximum no. of iterations allowed
tol   = 1e-8;
maxit = 20;

% initialise waitbar
%wb=waitbar(0,'Calculation running...');

% 1.) given (scalar) e_(n+1) at time t_(n+1) = stretch(n+1), 
% internal variables (= sdv) and deformation gradient (= e) at time t_n
% initialise partition of e
ebar = partition(e);

% loop over time-steps
for n=1:steps-1
    
    % display waitbar (CAUTION: slows down calculation!!!)
    %waitbar(n/(steps-1));
    
    % display current time step
    %disp(['n = ', num2str(n)]);

    % initialise sbar in order to enter the while-loop and set no. of
    % iterations to zero
    sbar = ones(5,1);
    iter = 0;
    
    % loop over Newton-iterations
    while norm(sbar)>tol % 6.) check convergence
        
        iter = iter+1;
        
        if iter>maxit
            %close(wb);
            error(['No convergence after ', num2str(maxit), ' iterations'])
        else

            % 2.) total deformation
            eneu = [e11(n+1) ebar(3) ebar(5);
                    ebar(3)  ebar(1) ebar(4);
                    ebar(5)  ebar(4) ebar(2)];

            % 3.) constitutive law: algorithmic stresses and moduli
            if matflag == 1
                [sv,E,sdvup] = materialbox_hooke(eneu,sdv,mat,mflag);
            elseif matflag == 2
                [sv,E,sdvup] = materialbox_kauderer(eneu,sdv,mat,mflag);
            elseif matflag == 3
                de = eneu - e;
                [sv,E,sdvup] = materialbox_maxwell(eneu,de,s,sdv,dt,mat,mflag);
%                [sv,E,sdvup] = materialbox_maxwell_parallel(eneu,de,s,sdv,dt,mat,mflag);
            elseif matflag == 4
                [sv,E,sdvup] = materialbox_vonmises(eneu,sdv,mat, mflag);
            end;

            % 4.) partitioning
            sbar = partition(sv);
            Ebar = partition(E);

            % 5.) update of lateral deformations
            ebar = ebar - Ebar\sbar;

            % display convergence
            disp(['iteration =', num2str(iter)])
            disp(['|sbar| = ', num2str(norm(sbar))])

        end % if

    end % while
    
    % update converged stress and strain tensor for next loading step
    s = [sv(1) sv(4) sv(6);
         sv(4) sv(2) sv(5);
         sv(6) sv(5) sv(3)];
     
    e = [e11(n+1) ebar(3) ebar(5);
         ebar(3)  ebar(1) ebar(4);
         ebar(5)  ebar(4) ebar(2)];
    
    % update internal variables
    sdv = sdvup;
    
    % store quantities for post-processing
    s11(n+1)=s(1,1);
    ep11(n+1) = sdv.ep(1,1);

end % for

%close(wb)

%% postprocessing ---------------------------------------------------------

figure('units', 'pixels', 'position', [100 100 1280 800]);
hold on;

% plot loading/unloading path
subplot(1,3,1)
hc = plot(time,e11);
set(hc                          , ...
   'LineWidth'       , 0.5 );
xlabel('time')
ylabel('stretch')
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [0 0 0], ...
  'YColor'      , [0 0 0], ...
  'LineWidth'   , 1         );
pbaspect([1 1 1])

% plot stress-stretch response
subplot(1,3,2)
hc = plot(e11,s11);
set(hc                          , ...
   'LineWidth'       , 0.5 );
xlabel('stretch')
ylabel('stress')
grid on
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [0 0 0], ...
  'YColor'      , [0 0 0], ...
  'LineWidth'   , 1         );
pbaspect([1 1 1])

% plot stress-stretch response
subplot(1,3,3)
hc = plot(time,s11);
set(hc                          , ...
   'LineWidth'       , 0.5 );
xlabel('time')
ylabel('stress')
grid on
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [0 0 0], ...
  'YColor'      , [0 0 0], ...
  'LineWidth'   , 1         );
pbaspect([1 1 1])
