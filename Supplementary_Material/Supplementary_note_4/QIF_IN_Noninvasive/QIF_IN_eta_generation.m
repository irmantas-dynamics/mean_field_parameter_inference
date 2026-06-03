function QIF_IN_eta_generation(N,fnm)
%
%   Generation of random eta from Lorentzian distribution with width delta
%   and peak at etab
%
%   Input: 
%           N   - number of neurons
%           fnm - number of iteration 
%
%   Output:
%           Script saves *.txt file with name pattern:
%           "eta_N%d_%d.mat", where the first %d is number of neurons and
%           second %d is number of iteration
%           This file contains time series of macroscopic variables:
%               Rt : spiking rate
%               v : mean potential
%               St: synaptic activation 
%           
%               PT: structure of neuron parameters
%               dt: integration step
%               tm: array of time
%               
%
    close all
    
    % model parameters
    etab = 4;
    eta=etab;
    delta = 0.3;

    
    % Very large values of eta near the distribution ends is ignored by
    % introducing small epss value
    epss = 1e-3;
    a_eta = eta+delta*tan(pi*(-0.5+epss));
    b_eta = eta+delta*tan(pi*(0.5-epss));
    % eta values for individual neurons    
    etam=truncated_cauchy_rand(eta, delta, a_eta, b_eta, N, 1);%(eta+delta*tan(pi*pag)).';

    fname = sprintf("time_series_rand/eta_N%d_%d.txt",N,fnm);
    writematrix(etam,fname);


end



function x = truncated_cauchy_rand(eta, delta, a, b, m, n)
    % eta   - location
    % delta - scale (>0)
    % [a,b] - truncation interval

    if delta <= 0
        error('delta must be positive');
    end
    if a >= b
        error('Require a < b');
    end

    % CDF of Cauchy
    F = @(x) (1/pi) * atan((x - eta)/delta) + 0.5;

    Fa = F(a);
    Fb = F(b);

    U = rand(m, n);
    U_scaled = Fa + (Fb - Fa) .* U;

    x = eta + delta * tan(pi * (U_scaled - 0.5));
end

