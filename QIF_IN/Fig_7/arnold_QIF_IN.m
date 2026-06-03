
  close all
  clear all

  dt= 0.01;

  
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

    
    dgkl = [0.95;0.95;0.95;0.95;0.95]*0+1;

    
     P.delta = dgkl(1)*PT.delta;
     P.J = dgkl(2) * PT.J;
     P.etab = dgkl(3) * PT.etab;
     P.taua = dgkl(4)*PT.taua;
     P.taum = dgkl(5)*PT.taum;

%     P = Pfin;

    Tosc= 27.56;
    KK.Amp = 2;
    KK.Tper = Tosc;
    KK.om = 2*pi/KK.Tper;
    
    Ttr= KK.Tper*100; % Transient time
    Dt = 3*KK.Tper; % gradient evaluation time
    Niter = 1000;
    Niter2 = Niter;
    
    tint  = linspace(0,round(5000/KK.Tper)*KK.Tper,50000);
    opts = odeset('RelTol',1e-7,'AbsTol',1e-9);

    Namp = 200;
    ampM = linspace(-1,1,Namp);
    NT = 200;
    TperM = linspace(0.7,1.3,NT)*Tosc;

    dist = NaN(Namp,NT);
    LM = NaN(Namp,NT);


    yinit = abs(rand(6,1));
    yinit0 =yinit;



   parfor i1 = 1:Namp
        
        K =KK;
        K.Amp = ampM(i1);
        K.Tper = TperM(1);
        K.om = 2*pi/TperM(1);
        disttmp = NaN(1,NT);
        LMtmp = NaN(1,NT);
        yinit =yinit0;
        [Tc, Yc] = ode45(@(t,x)orig_mod(t,x,PT,P,K), [0,Ttr],[yinit;0],opts);
        yinit = Yc(end,1:6).';
        for i2 = 1:NT
            K.Tper = TperM(i2);
            K.om = 2*pi/TperM(i2);
            [Tc, Yc] = ode45(@(t,x)orig_mod(t,x,PT,P,K), [0,100*K.Tper],[yinit;0],opts);
            yinit = Yc(end,1:6).';
    
            [Tc, Yc] = ode45(@(t,x)orig_mod(t,x,PT,P,K), [0,K.Tper],[yinit;0],opts);
            yinit = Yc(end,1:6).';
        
            iniVec = Yc(1,1:3);% - Yc(1,4:6);
            finVec = Yc(end,1:3);% - Yc(end,4:6);
            
            disttmp(1,i2) = norm(iniVec-finVec);
            LMtmp(1,i2) = Yc(end,end)/Tc(end);
        end
        dist(i1,:) = disttmp;
        LM(i1,:) = LMtmp;

   end



    figure
    imagesc(TperM,ampM,dist)
    title('norm[X(T)-X(0)]')
    xlabel('T_{ext}')
    ylabel('Amp')
    set(gca,'YDir','normal')
    clim([0 0.05])

    save('arnold_QIF_IN.mat','TperM','ampM','dist','PT','Tosc')
    

          


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


function dy =orig_mod(t,x,PT,P,K)


    v = x(2);

    vt = x(5);
    
    dy=zeros(4,1);

    Iext = K.Amp*(1+sin(K.om*t)/2)^3;
%     Iext = K.Amp*sin(K.om*t);

    % MODEL
    dy(1:3) = eqsns(t,x(1:3),P);
    dy(2)= dy(2) + Iext/P.taum; % it is necessary to norm Iext to time scale tau if such exists;

    
    % TRUE SYSTEM
    dy(4:6) = eqsns(t,x(4:6),PT);
    dy(5)= dy(5) + Iext/PT.taum;


    % loss function
    dy(7) = (v - vt)^2;

end




