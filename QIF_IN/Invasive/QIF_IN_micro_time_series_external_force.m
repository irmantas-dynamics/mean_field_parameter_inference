function QIF_IN_micro_time_series_external_force(N)
%
% Modeling network of QIF neurons with inhibitorry interaction, 
% according to "Firing rate equations require a spike synchrony mechanism
% to correctly describe fast oscillations in inhibitory networks" by
% F. Devalle, A. Roxin and E. Montbrió 
% https://doi.org/10.1371/journal.pcbi.1005881
% 
% Additiona periodic forcing is added.
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
%           "QIF_IN_N%d_ext.mat", where %d is number of neurons
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


    % parameters of external force
    K.Amp = -0.45;
    K.Tper = 28.;
    K.om = 2*pi/K.Tper;
    
    % Very large values of eta near the distribution ends is ignored by
    % introducing small epss value
    epss = 1e-3;
    pag=linspace(-0.5+epss,0.5-epss,N);
    etam=(etab+delta*tan(pi*pag)).';

    %initial values of the macroscopic variables
    V0 = 0;
    r0 = 0.1;

    j1=1:N;
    %initial values of potential of every neuron
    v0 =  (V0+r0*pi*tan(pi/2*(2*j1-N-1)/(N+1)))';


    %potential converted to theta variable
    tet = 2*atan(v0);
     st = 0.005;  


    P.N    =   N; 
    P.J    =   J; 
    P.etam = etam;
    P.taum = taum;
    P.taua = taua;
    P.etab = etab;
    P.delta = delta;

    % integration time
    Tend= 100*K.Tper;
    dt = 1e-3*taum;
    Nt = round(Tend/dt)+1;

    Rt = zeros(Nt,1);

    tm=linspace(0,(Nt-1)*dt,Nt);    
    Zm=zeros(Nt,1);
    St= zeros(Nt,1);

    y= [tet;st];


    tic
    t=0;
    fprintf('Initial integration\n')
    for i1=1:Nt
        
        if mod(i1/(Nt-1)*100,10)==0
            fprintf("%.1f proc.\n",i1/(Nt-1)*100)

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
        St(i1) = mean(y(N+1));
        Rt(i1) =  1/(pi*taum)*real((1-conj(Zm(i1)))./(1+conj(Zm(i1))));


        y =rk4(@(t,y)eqsns_micro(t,y,P,K),t,y,dt);
        y(1:N) = mod(y(1:N),2*pi);
        t = t+dt;
    end
    toc

    fprintf('End of microscopic equations integration\n\n\n')

    % mean potential and spiking rates estimated from Kuramoto order
    % parameter
    v =imag((1-conj(Zm))./(1+conj(Zm)));

    
    figure(2)
    subplot(311), hold on
    plot(tm,Rt)
    ylabel('r')
    subplot(312), hold on
    plot(tm,v), hold on
    ylabel('v')
    subplot(313), hold on
    plot(tm,St), hold on
    ylabel('S')

    PT.etab = etab;
    PT.delta = delta;
    PT.J = J;
    PT.taum =taum;
    PT.taua = taua;

    pav = sprintf('time_series/QIF_IN_N%d_ext.mat',N);
    
    save(pav,'PT','v','tm','tet','St','Nt','Rt','dt')   

%     st = y(N+1);
%     tet = y(1:N);
%     save('init.mat','tet','st')

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
% inhibitory interaction
        N    = P.N;
        J    = P.J;
        etam = P.etam;
        taum = P.taum;
        taua = P.taua;

        Zm = mean(exp(1i*y(1:N)));
        Zb= conj(Zm);
        W = (1-Zb)./(1+Zb);   
        R  = 1/(pi*taum)*real(W);

        Iext = K.Amp*(1+sin(K.om*t)/2)^3;
        
        cs=cos(y(1:N));
        st = y(N+1);


        dy = zeros(N+1,1);
        dy(1:N)     = (  1-cs+(1+cs).*(etam - J*taum*st +Iext ) )/taum;
        dy(N+1) = ( -y(N+1) +  R )/taua;
end
