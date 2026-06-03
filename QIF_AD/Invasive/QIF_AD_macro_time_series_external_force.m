function QIF_AD_macro_time_series_external_force()
% Modeling macroscopic equations of QIF neuron network with spike frequency
% adaptation, according to "Low-dimensional model for adaptive
% networks of spiking neurons" by  B. Pietras, P. Clusella and E. Montbrió 
% https://doi.org/10.1103/PhysRevE.111.014422
% Additional periodic external force is added
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



    K.Amp = -4;
    K.Tper = 80;
    K.om = 2*pi/K.Tper;



    % transient integration time
    T0= 60*K.Tper;
    % integration time for signal recording
    Tend= 100*K.Tper;
    %integration step
    dt = 1e-2;

    %number of integration steps
    Nt = round(Tend/dt)+1;

    tm = linspace(0,(Nt-1)*dt,Nt);

    yinit = [0.1;0.5;0.3];

    options = odeset('AbsTol',1e-9,'RelTol',1e-9);

    [T, Y] = ode45(@(t,x)orig_mod(t,x,P,K), [0,T0],yinit,options);
    yinit = Y(end,1:3).';

    [T, Y] = ode45(@(t,x)orig_mod(t,x,P,K), tm,yinit,options);


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
    pav = sprintf('time_series/QIF_AD_Ninf_chaos_ext.mat');
    save(pav,'tm','v','r','Aj','PT','K');
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
    
    dy=zeros(3,1);

    Iext = K.Amp*(1+sin(K.om*t)/2)^3;
%     Iext = K.Amp*sin(K.om*t);

    % MODEL
    dy(1:3) = eqsns(t,x(1:3),P);
    dy(2)= dy(2) + Iext/taum; % it is necessary to norm Iext to time scale tau if such exists;
    dy(3)= dy(3) + Iext/taua; % it is necessary to norm Iext to time scale tau if such exists;

end




