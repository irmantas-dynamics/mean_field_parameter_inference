function lyapunov_QIF_IN_K()
close all

    nit = 12*6;
    KM= linspace(0,2,nit);
    LPM = NaN(3,nit);
    maxL = NaN(1,nit);
    

    parfor i1=1:length(KM)
        LPM(:,i1)=lyapunov_QIF_IN(KM(i1));
        maxL(i1) = max( LPM(:,i1));
    end


    save('lyap_QIF_IN.mat','maxL','LPM','KM')

    figure
    plot(KM,LPM(1,:)), hold on
    plot(KM,LPM(2,:))
    plot(KM,LPM(3,:))
    plot(KM,KM*0,':')
    xlabel('K')
    ylabel('Lyap exp')

    figure
    plot(KM,maxL)



end

function LP=lyapunov_QIF_IN(K)

    % Lyapunov exponent calculation for the Montbrio system
    % Standard Gram-Schmidt orthogonalization algorithm

    % Model parameters
    etab = 4;       % Baseline excitability parameter
    delta = 0.3;    % Half-width of Lorentzian distribution of neuronal thresholds
    J = 21;         % Synaptic coupling strength
    taua = 5;       % Synaptic time constant
    taum = 10;      % Membrane time constant


    % Pack parameters into struct for easy passing to ODE functions
    P.taum  = taum;
    P.etab  = etab;
    P.J     = J;
    P.delta = delta;
    P.taua  = taua;

    % --- Find a point on (or near) the attractor via initial integration ---
    yinit = [0.1, 0, 0];                                    % Initial guess: [r, v, A]
    [T, Y] = ode45(@(t,y)eqsn(t,y,P), [0 5000], yinit);   % Integrate nonlinear system to t=5000
    yn = Y(end,:).';                                         % Use final state as starting point on attractor

    % --- Lyapunov integration settings ---
    tint  = 1000;           % Total integration time for Lyapunov computation
    stept = 0.051;          % Time interval between successive orthonormalizations
    nit   = round(tint/stept);  % Number of integration steps
    npral = 2000;           % Number of transient steps to skip before accumulating exponents

    % Initial conditions for the linearized (variational) equations:
    % 3 perturbation vectors stacked as a flattened 9-element column (identity matrix = orthonormal basis)
    yl = [1 0 0  0 1 0  0 0 1]';

    % Combined initial state: nonlinear state (3) + variational vectors (9)
    yinit = [yn; yl];

    fcn_integrator = @ode45;    % Choose ODE integrator

    % --- Transient (warm-up) phase: advance solution without accumulating exponents ---
    t = 0;
    for j = 1:npral
        % Integrate the extended system (nonlinear + 3 linearized copies) over one step
        sol = fcn_integrator(@(t,y)eqsn_ext(t,y,P,K), [t t+stept], yinit);
        Y = deval(t+stept, sol);    % Extract state at end of step
        t = t + stept;

        % Re-orthonormalize the perturbation vectors via Gram-Schmidt
        [x0, ~] = gram_sch(Y(4:12));   % x0 = orthonormal basis; norms discarded during warm-up

        % Update initial conditions for next step
        yinit(1:3)  = Y(1:3);   % Nonlinear state
        yinit(4:12) = x0;       % Orthonormalized perturbation vectors
    end

    % --- Main Lyapunov exponent computation loop ---
    t   = 0;
    lp  = zeros(3, nit);    % Matrix to store running Lyapunov exponent estimates [3 x nit]
    T   = zeros(1, nit);    % Time vector
    sum = zeros(3, 1);      % Cumulative sum of log(norms) for each exponent

    for j = 1:nit
        % Integrate extended system over one step
        sol = fcn_integrator(@(t,y)eqsn_ext(t,y,P,K), [t t+stept], yinit);
        Y = deval(t+stept, sol);
        t = t + stept;
        T(j) = t;

        % Gram-Schmidt: get orthonormal vectors x0 and their norms d
        [x0, d] = gram_sch(Y(4:12));

        % Accumulate log-norms: each norm reflects local stretching/contraction
        sum = sum + log(d);

        % Running average of Lyapunov exponents (sum of log-norms divided by elapsed time)
        lp(:, j) = sum / t;

        % Update initial conditions
        yinit(1:3)  = Y(1:3);
        yinit(4:12) = x0;
    end

    LP = lp(:,end);

%     % Display final Lyapunov exponent estimates
%     disp(lp(:, end))
% 
%     % Plot convergence of all three Lyapunov exponents over time
%     figure
%     plot(T, lp)
%     xlabel('Time'); ylabel('\lambda'); title('Lyapunov Exponents (Montbrio System)');
end


function [x0, d] = gram_sch(x)
    % Gram-Schmidt orthonormalization
    % Input:  x  - 9-element column vector (three stacked 3D perturbation vectors)
    % Output: x0 - 9-element orthonormalized column vector (same structure)
    %         d  - 3-element vector of norms before normalization

    d  = zeros(3, 1);   % Array to store vector norms

    x1 = x(1:3);        % First perturbation vector
    x2 = x(4:6);        % Second perturbation vector
    x3 = x(7:9);        % Third perturbation vector

    % Normalize first vector
    d(1)  = norm(x1);
    x10   = x1 / d(1);

    % Remove projection of x2 onto x10, then normalize
    a2    = x2 - (x2' * x10) * x10;
    d(2)  = norm(a2);
    x20   = a2 / d(2);

    % Remove projections of x3 onto x10 and x20, then normalize
    a3    = x3 - (x3' * x10) * x10 - (x3' * x20) * x20;
    d(3)  = norm(a3);
    x30   = a3 / d(3);

    % Stack orthonormal vectors back into a 9-element column
    x0 = [x10; x20; x30];
end


function dy = eqsn(t, y, P, K)
    % State variables: y(1) = r (mean firing rate)
    %                  y(2) = v (mean membrane potential)
    %                  y(3) = A (synaptic activation / adaptation variable)

    taum  = P.taum;
    etab  = P.etab;
    J     = P.J;
    delta = P.delta;
    taua  = P.taua;

    dy = zeros(3, 1);

    % Firing rate equation (r-dot): driven by delta and the r*v product
    dy(1) = (delta / (pi * taum) + 2 * y(1) * y(2)) / taum;

    % Mean potential equation (v-dot): includes quadratic nonlinearity and synaptic input
    dy(2) = (y(2)^2 - (taum * pi * y(1))^2 + etab - J * taum * y(3)) / taum;

    % Synaptic/adaptation variable equation (A-dot): first-order low-pass filter of r
    dy(3) = (-y(3) + y(1)) / taua;
end


function f = eqsn_ext(t, X, P, K)
    % Extended ODE system: nonlinear  system + three linearized (variational) copies
    % Used to simultaneously evolve the trajectory and the perturbation vectors
    % X(1:3)   = nonlinear state [r, v, A]
    % X(4:6)   = first perturbation vector
    % X(7:9)   = second perturbation vector
    % X(10:12) = third perturbation vector

    taum = P.taum;
    etab = P.etab;
    J    = P.J;
    taua = P.taua;

    % Unpack nonlinear state
    r = X(1); v = X(2); A = X(3);

    % Arrange three perturbation vectors as columns of a 3x3 matrix
    Y = [X(4:6), X(7:9), X(10:12)];

    f = zeros(12, 1);

    % Nonlinear system equations
    f(1:3) = eqsn(t, X(1:3), P);

    % Jacobian of the Montbrio system evaluated at the current state (r, v, A)
    % Rows correspond to d(r-dot)/d*, d(v-dot)/d*, d(A-dot)/d*
    Jac = [      2*v/taum,          2*r/taum,        0;
          -2*pi^2*taum*r,    2*v/taum - K,           -J;
                1/taua,              0,          -1/taua];

    % Linearized (variational) equations: dY/dt = Jac * Y
    % Flattened back to a 9-element column for the ODE integrator
    f(4:12) = Jac * Y;
end
