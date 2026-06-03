function QIF_AD_invasive_spectra_comparisson()
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



    K.Amp = 0;%-4;
    K.Tper = 80;
    K.om = 2*pi/K.Tper;

    NnM = {'200','500','1000','10000','100000','inf'};

%     NnM = {'100000'};

    plt = [];
    lgnd = {};

    clrs =  [0, 0.4470, 0.7410;...
            0.8500, 0.3250, 0.0980;...
            0.9290, 0.6940, 0.1250;...
    	    0.4940, 0.1840, 0.5560;...	
            0.4660, 0.6740, 0.1880;...
          	0.3010, 0.7450, 0.9330;...
            0.6350, 0.0780, 0.1840];

    for i1 = 1:length(NnM)
        Nn = NnM{i1};%'200'
    
        fname = sprintf("./time_series/QIF_AD_N%s_chaos.mat",Nn);
        load(fname);

        
        dt = tm(2)-tm(1);

    
    
        fs = 1/dt;
        x= v;
    
        nHann = 10;
    
        figure(1)
        subplot(3,2,i1)
        [Power, Frequency] = WelchPowerSpectralDensity(x, [], Hann(round(length(x) /nHann)), 0., fs); 
    %     plot(gca, Frequency, 10*log10(abs(Power)))
        plt(i1) = semilogy(gca, Frequency, (abs(Power)), Color=clrs(1,:));
        hold on
    
        fname2 = sprintf("./optimization_results/QIF_AD_%s_optim_results_ext.csv",Nn);
        Prm = readmatrix(fname2);
    
        
%         [~,idx] = min(Prm(:,end));
%         pp = Prm(idx,:);

        pp=median(Prm);
        

        P.delta = pp(1);
        P.J    = pp(2);
        P.etab = pp(3);
        P.beta = pp(4);
        P.taum = pp(5);
    


        pav = sprintf("reconstructed_model_time_series/signal_N%s.mat",Nn);
        load(pav,"Y","T");


        fs = 1/(T(2)-T(1));

    
    
        x= Y(:,2);
         nHann = 10;
    
         figure(1)
        [Power, Frequency] = WelchPowerSpectralDensity(x, [], Hann(round(length(x) /nHann)), 0., fs); 
        semilogy(gca, Frequency, (abs(Power)),Color=clrs(2,:),LineStyle='-',LineWidth=1.5)
        %semilogy(gca, Frequency, (abs(Power)))
        lgnd(2*i1-1) = {sprintf("N =%s",NnM{i1})};
        lgnd(2*i1) = {sprintf("%s - %s",NnM{i1},'pred')};
        xlim([0,2.6])
        ylim([1e-6 1e2])
        legend(plt(i1),lgnd(2*i1-1))
%         legend({NnM{i1},'fit'})
%         ylim([1e-7 1e3])
        if mod(i1,2)==1
            ylabel('Spectral density')
            yticks([1e-6 1e-4 1e-2 1e0 1e2])
        else
            yticks([1e-6 1e-4 1e-2 1e0 1e2])
            set(gca, 'YTickLabel', []);
        end

        if i1<5
             set(gca, 'XTickLabel', []);
        end

        if i1>4
            xlabel('Frequency, kHz')
        end
       


    end



end

