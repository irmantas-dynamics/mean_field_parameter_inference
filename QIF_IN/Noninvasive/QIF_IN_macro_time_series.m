function QIF_IN_macro_time_series()
%
% Modeling mean field equations for the thermodynamic limit of the
% network of QIF neurons with inhibitorry interaction, 
% according to "Firing rate equations require a spike synchrony mechanism
% to correctly describe fast oscillations in inhibitory networks" by
% F. Devalle, A. Roxin and E. Montbrió 
% https://doi.org/10.1371/journal.pcbi.1005881
%
%
% Integration is performed with Runge-Kutta algorithm
%
%   Input: 
%   Output:
%           Script saves *.mat file with name pattern:
%           "time_series/QIF_IN_Ninf.mat"
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

    % external force amplitude
    K=0;


    PT.taum = taum;
    PT.etab = etab;
    PT.J    = J   ;
    PT.delta = delta;
    PT.taua = taua;

    P = PT
    
    % period of the macroscopic variables
    Tosc = 27.71;

    % integration time
    Tend= 70*Tosc;
    %integration step
    dt = 1e-3*taum;t
    Nt = round(Tend/dt)+1;
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



    save('time_series/QIF_IN_Ninf.mat','tm','v','r','St','PT');

end


function y = rk4(F_xy,t,y,h)   
    k_1 = F_xy(t,y);
    k_2 = F_xy(t+0.5*h,y+0.5*h*k_1);
    k_3 = F_xy((t+0.5*h),(y+0.5*h*k_2));
    k_4 = F_xy((t+h),(y+k_3*h));
    y = y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
end

function dy = eqsn(t,y,P,K)

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

function dy = orig_mod(t,y,P,K,PT)

    taum = P.taum;

    dy=zeros(7,1);

    dy(1:3) = eqsn(t,y(1:3),P,K);
    dy(2) = dy(2) + K*(y(5)-y(2))/taum;
    dy(4:6) = eqsn(t,y(4:6),PT,K);
    dy(7) = 0.5*(y(2)-y(5))^2;

end
