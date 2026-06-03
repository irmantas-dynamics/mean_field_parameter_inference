function QIF_AD_micro_time_series(N)
%
% Modeling network of QIF neurons with spike frequency adaptation, according to "Low-dimensional model for adaptive
% networks of spiking neurons" by  B. Pietras, P. Clusella and E. Montbrió 
% https://doi.org/10.1103/PhysRevE.111.014422
%
% Equations of QIF neurons are converted to Theta neuron model.
% Mean potential and spiking rate is computed from order parameter using
% conformal mapping.
%
% Integration is performed with Runge-Kutta algorithm
%
%   Input: 
%           N - number of neurons
%
%   Output:
%           Script saves *.mat file with name pattern:
%           "QIF_AD_N%d_chaos.mat", where %d is number of neurons
%           This file contains time series of macroscopic variables:
%               r : spiking rate
%               v : mean potential
%               Aj: average adaptation variable
%           
%               PT: structure of neuron parameters
%               dt: integration step
%               tm: array of time
%               
%
    close all

    %amplitude of external force
    K= 0;

    % Neuron parameters
    eta = 3.25; % center of Lorenz distribution
    delta = 1;
    J = 20;
    taua = 100;
    taum = 10;
    beta = 1;

    
    % Very large values of eta near the distribution ends is ignored by
    % introducing small epss value
    epss = 1e-3;
    pag=linspace(-0.5+epss,0.5-epss,N);
    % eta values for individual neurons    
    etam=(eta+delta*tan(pi*pag)).';

    %initial values of the macroscopic variables
    V0 = 0;
    r0 = 0.1;

    j1=1:N;
    %initial potential of every neuron
    v0 = (V0+r0*pi*taum*tan(pi/2*(2*j1-N-1)/(N+1)))'; 

    %potential converted to theta variable
    tet = 2*atan(v0);
    aj = zeros(N,1);

    % load('init.mat','tet','aj')

    P.N    =   N; 
    P.J    =   J; 
    P.etam = etam;
    P.beta = beta;
    P.taum = taum;
    P.taua = taua;
    P.etab = eta;
    P.delta = delta;

    %integration time of transitional dynamics
    T0= 1000;
    %integration time for signal recording
    Tend= 5000;

    %integration step
    dt = 1e-2;

    %number of integration steps
    Nt = round(Tend/dt)+1;
    Nt0 = round(T0/dt)+1;
    

    tm = linspace(0,(Nt-1)*dt,Nt);
    %arrays for storing Kuramoto order parameter 
    Zm = zeros(Nt,1);
    % mean adaptaion variable
    Aj = zeros(Nt,1);

       
    %load('init.mat','tet','aj')
    %initial condition for the microscopic equations
    y = [tet;aj];


    t = 0;
    tic
    fprintf('Initial integration\n')
    for i1=1:Nt0
        
        if mod(i1/(Nt0-1)*100,10)==0
            fprintf("%.1f proc.\n",i1/(Nt0-1)*100)

        end
        y =rk4(@(t,y)eqsns_micro(t,y,P,K),t,y,dt);
        y(1:N) = mod(y(1:N),2*pi);
        t = t+dt;
    end
    toc


    t = 0;
    tic
    fprintf('\n\n\nSignal generation\n')
    for i1=1:Nt
        
        if mod(i1/(Nt-1)*100,10)==0
            fprintf("%.1f proc.\n",i1/(Nt-1)*100)

        end

        Zm(i1) = mean(exp(1i*y(1:N)));
        Aj(i1) = mean(y(N+1:2*N));

        y =rk4(@(t,y)eqsns_micro(t,y,P,K),t,y,dt);
        y(1:N) = mod(y(1:N),2*pi);
        t = t+dt;
    end
    toc

    fprintf('End of microscopic equations integration\n\n\n')

    % mean potential and spiking rates estimated from Kuramoto order
    % parameter
    v =imag((1-conj(Zm))./(1+conj(Zm)));
    r =1/(pi*taum)*real((1-conj(Zm))./(1+conj(Zm)));

   
    
    figure(1)
    subplot(311)
    plot(tm,r), hold on
    ylabel('r')
    subplot(312)
    plot(tm,v), hold on
    ylabel('v')
    subplot(313)
    plot(tm,Aj), hold on
    ylabel('A')

    PT = P;

    pav = sprintf('time_series/QIF_AD_N%d_chaos.mat',N);
    
    save(pav,'PT','v','r','tm','tet','Aj','Nt','dt')

    %aj = y(N+1:2*N);
    %tet = y(1:N);
    %save('init.mat','tet','aj')

end


    

function y = rk4(F_xy,t,y,h)   
    k_1 = F_xy(t,y);
    k_2 = F_xy(t+0.5*h,y+0.5*h*k_1);
    k_3 = F_xy((t+0.5*h),(y+0.5*h*k_2));
    k_4 = F_xy((t+h),(y+k_3*h));
    y = y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  
end

function dy = eqsns_micro(t,y,P,K)  
% Microscopic network simulations of QIF neurons with
% quadratic spike-frequency adaptation 
        N    = P.N;
        J    = P.J;
        etam = P.etam;
        beta = P.beta;
        taum = P.taum;
        taua = P.taua;

        Zm = mean(exp(1i*y(1:N)));
        Zb= conj(Zm);
        W = (1-Zb)./(1+Zb);   
        R  = 1/(pi*taum)*real(W);

        Iext = 0;%K.Amp*(1+sin(K.om*t)/2)^3;
        
        cs=cos(y(1:N));
        Fn = (etam - y(N+1:2*N) + J*taum*R+Iext);

        dy = zeros(2*N,1);
        dy(1:N)     = ( 1-cs+(1+cs).*Fn )/taum;
        dy(N+1:2*N) = ( -y(N+1:2*N) +   beta*Fn )/taua;
end



