close all
clear all

dt= 0.01;

etab = 3.2;
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


PT =P;
P.Jdaug = 10;
PT.J = PT.J/P.Jdaug;
PT.Jdaug= P.Jdaug;


KK.Amp = 2;
KK.om = 2*pi/100;
KK.Tper = 2*pi/KK.om;

Ttr= KK.Tper*100; % Transient time
Dt = 3*KK.Tper; % gradient evaluation time
Niter = 1000;
Niter2 = Niter;

tint  = linspace(0,round(5000/KK.Tper)*KK.Tper,50000);
opts = odeset('RelTol',1e-7,'AbsTol',1e-9);

Namp = 200;
ampM = linspace(-5,5,Namp);
NT = 200;
TperM = linspace(20,110,NT);



dist = NaN(Namp,NT);
LM = NaN(Namp,NT);



yinit0 = randn(3,1)+1;


yinit =yinit0;
[Tc, Yc] = ode45(@(t,x)eqsns(t,x,P,KK), [0,Ttr],yinit,opts);



tic
parfor i1 = 1:Namp
    K =KK;
    K.Amp = ampM(i1);
    K.Tper = TperM(1);
    K.om = 2*pi/TperM(1);
    disttmp = NaN(1,NT);
    yinit =yinit0;
    [Tc, Yc] = ode45(@(t,x)eqsns(t,x,PT,K), [0,Ttr],yinit,opts);
    yinit = Yc(end,1:3).';
    for i2 = 1:NT
        K.Tper = TperM(i2);
        K.om = 2*pi/TperM(i2);
        [Tc, Yc] = ode45(@(t,x)eqsns(t,x,PT,K), [0,50*K.Tper],yinit,opts);
        yinit = Yc(end,1:3).';

        [Tc, Yc] = ode45(@(t,x)eqsns(t,x,PT,K), [0,K.Tper],yinit,opts);
        yinit = Yc(end,1:3).';
    
        iniVec = Yc(1,1:3);
        finVec = Yc(end,1:3);
       
        
        disttmp(1,i2) = norm(iniVec-finVec);

    end
    dist(i1,:) = disttmp;
end
toc

figure
imagesc(TperM,ampM,dist)
title('norm[X(T)-X(0)]')
xlabel('T_{force}')
ylabel('Amp')
set(gca,'YDir','normal')
%     clim([0 0.05])

save('arnold_QIF_AD.mat','TperM','ampM','dist','PT')



function dy = eqsns(t,y,P,K)

    taum = P.taum;
    etab = P.etab;
    J = P.J;
    Jdaug = P.Jdaug;
    delta = P.delta;
    beta = P.beta;
    taua = P.taua;

    dy = zeros(3,1);
    
%     Iext = K.Amp*sin(K.om*t);
    Iext = K.Amp*(1+sin(K.om*t)/2)^3;

    %spiking rate R
    dy(1) = ( delta/((1+beta)*pi*taum)+2*y(1)*y(2))/taum;
    %potencial
    dy(2) = ( y(2)^2 - (taum*pi*y(1))^2 + etab + Iext + J*Jdaug *taum*y(1)-y(3) )/taum;
    % adaptation variable
    dy(3) = (-y(3)*(1+beta)+beta*(etab+J*Jdaug*taum*y(1))+Iext)/taua;

end






