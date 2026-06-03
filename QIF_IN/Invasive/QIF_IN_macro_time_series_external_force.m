function QIF_IN_macro_time_series_external_force()
%
% Modeling network of QIF neurons with inhibitorry interaction, 
% according to "Firing rate equations require a spike synchrony mechanism
% to correctly describe fast oscillations in inhibitory networks" by
% F. Devalle, A. Roxin and E. Montbrió 
% https://doi.org/10.1371/journal.pcbi.1005881
%
% Equations of QIF neurons are converted to Theta neuron model.
% Mean potential and spiking rate is computed from order parameter using
% conformal mapping.
%
% Integration is performed with Runge-Kutta algorithm
%
%   Input: 
%
%   Output:
%           Script saves *.mat file with name pattern:
%           "QIF_IN_Ninf_ext.mat",
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
    delta = 0.3;
    J = 21;
    taua = 5;
    taum = 10;

    PT.taum = taum;
    PT.etab = etab;
    PT.J    = J   ;
    PT.delta = delta;
    PT.taua = taua;

    P = PT;

    % external force parameters
    ampM = -0.45;
    TperM = 28.;        

    K.Amp = ampM;
    K.Tper = TperM;
    K.om = 2*pi/TperM;

    % integration time
    Tend= 100*K.Tper;
    %integration step
    dt = 1e-3*taum;
    Nt = round(Tend/dt);
    tm=linspace(0,(Nt-1)*dt,Nt); 

    r = zeros(Nt,1);
    v = zeros(Nt,1);
    % mean adaptaion variable
    St = zeros(Nt,1);
    
    y = abs(randn(1,6));

    t=0;
    fprintf('Initial integration\n')
    for i1=1:Nt
        
        if mod(i1/(Nt-1)*100,10)==0
            fprintf("%f proc.\n",i1/(Nt-1)*100)

        end
        y =rk4(@(t,y)orig_mod(t,y,P,K,PT),t,y,dt);
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

        r(i1) = y(1);
        v(i1) = y(2);
        St(i1) = y(3);

        y =rk4(@(t,y)orig_mod(t,y,P,K,PT),t,y,dt);
        t = t+dt;
    end
    toc

    figure
    subplot(3,1,1)
    plot(tm,r), hold on
    subplot(3,1,2)
    plot(tm,v), hold on
    subplot(3,1,3)
    plot(tm,St), hold on
    
    save('time_series/QIF_IN_Ninf_ext.mat','tm','v','r','St','PT');


end


function y = rk4(F_xy,t,y,h)   
    k_1 = F_xy(t,y);
    k_2 = F_xy(t+0.5*h,y+0.5*h*k_1);
    k_3 = F_xy((t+0.5*h),(y+0.5*h*k_2));
    k_4 = F_xy((t+h),(y+k_3*h));
    y = y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
end


function dy = eqsns(t,y,P)

    taum = P.taum;
    etab = P.etab;
    J = P.J;
    delta = P.delta;

    taua = P.taua;

    dy = zeros(3,1);

    %spiking rate R
    dy(1) = ( delta/(pi*taum)+2*y(1)*y(2))/taum;
    %potencialas
    dy(2) = ( y(2)^2 - (taum*pi*y(1))^2 + etab - J *taum*y(3) )/taum;
    % synapse variable
    dy(3) = ( -y(3)+y(1) )/taua;

end


function dy =orig_mod(t,x,PT,K,P)
    
    dy=zeros(3,1);

    Iext = K.Amp*(1+sin(K.om*t)/2)^3;
%     Iext = K.Amp*sin(K.om*t);

    % MODEL
    dy(1:3) = eqsns(t,x(1:3),P);
    dy(2)= dy(2) + Iext/P.taum; % it is necessary to norm Iext to time scale tau if such exists;

end




