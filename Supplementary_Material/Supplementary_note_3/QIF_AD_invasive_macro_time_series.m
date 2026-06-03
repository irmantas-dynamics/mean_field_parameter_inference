function QIF_AD_invasive_macro_time_series()
% Modeling macroscopic equations of QIF neuron network with spike frequency
% adaptation, according to "Low-dimensional model for adaptive
% networks of spiking neurons" by  B. Pietras, P. Clusella and E. Montbrió 
% https://doi.org/10.1103/PhysRevE.111.014422
%
% System parameters are taken from optimization with invasive force results
% in folder 'optimization_results'.
%
%   Output:
%           Script saves *.mat file with name pattern:
%               "treconstructed_model_time_series/signal_N%d.mat"
%           This file contains time series of macroscopic variables:
%               r : spiking rate
%               v : mean potential
%               Aj: average adaptation variable
%           
%               PT: structure of neuron parameters
%               tm: array of time
%               
%

    close all

    etab = 3.25;
    delta= 1;    
    J    = 20;
    taua = 100;
    taum = 10;
    beta = 1;    


    P.taum = taum;
    P.etab = etab;
    P.J    = J   ;
    P.delta = delta;
    P.beta = beta;
    P.taua = taua;


    K.Amp = 0;
    K.Tper = 80;
    K.om = 2*pi/K.Tper;

    NnM = {'200','500','1000','10000','100000','inf'};


    parfor i1 = 1:length(NnM)


        Nn = NnM{i1};
        PP =P;

        % transient integration time
        T0= 60*K.Tper;
        % integration time for signal recording
        Tend= 500;%10*K.Tper;
        %integration step
        dt = 1e-2;
    

        fname2 = sprintf("./optimization_results/QIF_AD_%s_optim_results_ext.csv",Nn);
        Prm = readmatrix(fname2);
     

        pp=median(Prm);
        

        PP.delta = pp(1);
        PP.J    = pp(2);
        PP.etab = pp(3);
        PP.beta = pp(4);
        PP.taum = pp(5);
    
        yinit = rand(3,1);
        [T, Y] = eulermaruyama(@(t,x)orig_mod(t,x,PP,K),@noise, 0:dt:T0,yinit,0,dt*1e-2);
        yinit = Y(end,1:3).';
    
        [T, Y] = eulermaruyama(@(t,x)orig_mod(t,x,PP,K),@noise, 0:dt:Tend,yinit,0,dt*1e-2);


        pav = sprintf("reconstructed_model_time_series/signal_N%s.mat",Nn);
        s =struct("T",T,"Y",Y);
        save(pav,"-fromstruct",s);

    end




end

function dg = noise(t,y)
    dg = [1,1,1]*1e-6;%[0.00025,0.000,0];
end

function dy = eqsns(t,y,P)

    taum = P.taum;
    etab = P.etab;
    J = P.J;
    delta = P.delta;
    beta = P.beta;
    taua = P.taua;

    dy = zeros(3,1);



    %spiking rate R
    dy(1) = ( delta/((1+beta)*pi*taum)+2*y(1)*y(2))/taum;
    %potencialas
    dy(2) = ( y(2)^2 - (taum*pi*y(1))^2 + etab + J *taum*y(1)-y(3) )/taum;
    % adaptation variable
    dy(3) = (-y(3)*(1+beta)+beta*(etab+J*taum*y(1)))/taua;

end


function dy =orig_mod(t,x,P,K)

    taum = P.taum;
    taua = P.taua;
    
    dy=zeros(1,3);

    Iext = 0;
%     Iext = K.Amp*sin(K.om*t);

    % MODEL
    dy(1:3) = eqsns(t,x(1:3),P);
    dy(2)= dy(2) + Iext/taum; % it is necessary to norm Iext to time scale tau if such exists;
    dy(3)= dy(3) + Iext/taua; % it is necessary to norm Iext to time scale tau if such exists;

end

function [tspan,x] = eulermaruyama(f,L,tspan,x0,Q,dt)
%% EULERMARUYAMA - Numerical SDE solver: The Euler-Maruyama method
%
% Syntax:
%   [tspan,x] = eulermaruyama(f,L,tspan,x0,Q,dt)
%
% In:
%   f      - Drift function, f(t,x)
%   L      - Diffusion function, L(t,x)
%   tspan  - Output time steps, [t0,...,tend]
%   x0     - Initial condition
%   Q      - Spectral density (default: standard Brownian motion)
%   dt     - Fine integration step size (default: tspan(2)-tspan(1))
%
% Out:
%   tspan  - Output time steps
%   x      - Solved values at tspan steps
%
% Description:
%   Integrates the system of stochastic differential equations
%     dx = f(x,t) dt + L(x,t) dbeta,  for x(0) = x0
%   over the time interval defined in tspan. The output is returned at the
%   coarse tspan steps, but the internal integration is performed using the
%   finer step size dt.
%
% Copyright:
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying
%   LICENSE file for details.

%%

  % Check if Q given
  if nargin < 5 || isempty(Q), Q = eye(size(L(tspan(1),x0),2)); end

  % Default fine step equals the coarse step (original behaviour)
  DT = tspan(2) - tspan(1);
  if nargin < 6 || isempty(dt), dt = DT; end

  % Validate dt
  if dt > DT
    error('Fine integration step dt must not exceed the output step DT.');
  end

%   % Cholesky factor of Q
%   cQ = chol(Q,'lower');

  % Number of output steps
  steps = numel(tspan);

  % Allocate output
  x = zeros(steps, numel(x0));

  % Initial state
  x(1,:) = x0(:);

  % Iterate over output intervals
  for k = 2:steps

    % Coarse output interval
    T0 = tspan(k-1);
    T1 = tspan(k);

    % Build fine sub-grid for this interval.
    % Use as many full dt steps as fit, then one remainder step if needed.
    n_full = floor((T1 - T0) / dt);
    t_fine = T0 + (0:n_full) * dt;
    if t_fine(end) < T1               % add remainder node if needed
      t_fine(end+1) = T1; %#ok<AGROW>
    end

    % Integrate over the sub-grid
    xk = x(k-1,:);
    for j = 2:numel(t_fine)

      dts  = t_fine(j) - t_fine(j-1);
      db   = sqrt(dts) * randn(size(Q,1), 1);
      xk   = xk + ...
             f(t_fine(j-1), xk) * dts + ...
             L(t_fine(j-1), xk) * db;

    end

    % Store output at the coarse step
    x(k,:) = xk';

  end

end



