function QIF_AD_macro_time_series()
% Modeling macroscopic equations of QIF neuron network with spike frequency
% adaptation, according to "Low-dimensional model for adaptive
% networks of spiking neurons" by  B. Pietras, P. Clusella and E. Montbrió 
% https://doi.org/10.1103/PhysRevE.111.014422
%
%   Output:
%           Script saves *.mat file with name pattern:
%               "time_series/QIF_AD_Ninf_chaos.mat"
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

    % transient integration time
    T0= 1000;
    % integration time for signal recording
    Tend= 5000;

    %integration step
    dt = 1e-2;

    %number of integration steps
    Nt = round(Tend/dt)+1;

    tm = linspace(0,(Nt-1)*dt,Nt);

    yinit = [ 0.1,   0.,  0];

    options = odeset('AbsTol',1e-9,'RelTol',1e-9);

    [T,Y] = ode45(@(t,y)eqsns(t,y,P),[0 T0],yinit,options);
    yinit= Y(end,:);

    [T,Y] = ode45(@(t,y)eqsns(t,y,P),tm,yinit,options);


    figure(2)
    subplot(311), hold on
    plot(T,Y(:,1))
    title('R')
    subplot(312), hold on
    plot(T,Y(:,2))
    title('V')
    subplot(313), hold on
    plot(T,Y(:,3))
    title('A')

    PT =P;
    r= Y(:,1);
    v= Y(:,2);
    Aj = Y(:,3);
    pav = sprintf('time_series/QIF_AD_Ninf_chaos.mat');
    save(pav,'PT','v','r','tm','Aj','Nt','dt')
end

function dy = eqsns(t,y,P)

    taum = P.taum;
    etab = P.etab;
    J = P.J;
    delta = P.delta;
    beta = P.beta;
    taua = P.taua;

    dy = zeros(3,1);

    %potencialas
    dy(1) = ( delta/((1+beta)*pi*taum)+2*y(1)*y(2))/taum;
    %spiking rate R
    dy(2) = ( y(2)^2 - (taum*pi*y(1))^2 + etab + J *taum*y(1)-y(3) )/taum;
    % adaptation variable
    dy(3) = (-y(3)*(1+beta)+beta*(etab+J*taum*y(1)))/taua;

end