close all
clear all

% actual parameters
prm = [0.3;21; 4;10;5];

bxWidth = 0.2;
deviation = 0.15; 

    

pv = {'$\delta \Delta$','$\delta J$','$\delta \bar{\eta}$','$\delta \tau_m$','$\delta \tau_a$'};
lbl = {'200','500','1000','10000','100000','inf'};
lbl2 = {'2\cdot 10^2','5\cdot 10^2','10^3','10^4','10^5','\infty'};

mnm = ones(length(lbl),1)*10000;
mxm = ones(length(lbl),1)*-10000;


dirname = 'Invasive/optimization_results';


ax1 = [0.07, 0.72,   0.40,  0.26;
       0.56, 0.72,   0.4,   0.26;
       0.07, 0.44,   0.40,  0.26;
       0.56, 0.44,   0.4,   0.26;
       0.27, 0.10,   0.40,  0.26]

% loop through all files with optimization results with periodic force
for ifile = 1 :length(lbl)

    fname = sprintf(strcat(dirname,'/QIF_IN_%s_ext_optim_results.csv'),lbl{ifile})
    M = readmatrix(fname);
    tmp =ones(size(M,1),1);

    % Figure of parameters
    figure(1)
    NValues = NaN(size(M,1),length(lbl));
    % loop through parameters
    for i1 = 1:5
        
        ax2(i1)=subplot(3,2,i1);

        % parameter relative deviation
        tmpY = (M(:,i1)-prm(i1))/prm(i1);
        
        inBxChartPar(i1)=boxchart(tmp*ifile+ deviation,  tmpY);
        hold on
        inBxChartPar(i1).BoxFaceColor = [0.9 0.25 0.25];
        inBxChartPar(i1).BoxWidth = bxWidth;
        inBxChartPar(i1).MarkerStyle = '.';
        inBxChartPar(i1).MarkerColor = [0.9 0.25 0.25];
        ylabel(pv{i1},'Interpreter','latex')  
        
        tmpmnm = min(tmpY);
        tmpmax = max(tmpY);
        if mnm(i1)>tmpmnm
            mnm(i1) = tmpmnm;
        end
        if mxm(i1)<tmpmax
            mxm(i1) = tmpmax;
        end
           

    end

    set(gca,'XTickLabel',lbl)
%         xticks(1:6)

    % Figure of loss function
    figure(2)
    inBxChartLoss(ifile) = boxchart(tmp*ifile+deviation,  M(:,end)+1e-7); hold on
    inBxChartLoss(ifile).BoxFaceColor = [0.9 0.25 0.25];
    inBxChartLoss(ifile).BoxWidth = bxWidth;
    inBxChartLoss(ifile).MarkerStyle = '.';
    inBxChartLoss(ifile).MarkerColor = [0.9 0.25 0.25];
%     xticks(1:6)
    xlabel('$N$','Interpreter','latex')
    ylabel('Loss', 'Interpreter','latex')
    box on

end


  

mnm = ones(length(lbl),1)*10000;
mxm = ones(length(lbl),1)*-10000;

% fname = sprintf('theta_mean_%s_results.csv',lbl{1})
% M = readmatrix(fname);

dirname = 'Noninvasive/optimization_results';

% loop through all files with optimization results with feedback
for ifile = 1 :length(lbl)
    fname = sprintf(strcat(dirname,'/QIF_IN_%s_optim_results.csv'),lbl{ifile})
    M = readmatrix(fname);
    tmp =ones(size(M,1),1);
    
    % Figure of parameters
    figure(1)
    NValues = NaN(size(M,1),length(lbl));
    % loop through parameters
    for i1 = 1:5
        subplot(3,2,i1)

        % parameter relative deviation
        tmpY = (M(:,i1)-prm(i1))/prm(i1);

        fdBxChart(i1)=boxchart(tmp*ifile-deviation,  tmpY);
        hold on
        fdBxChart(i1).BoxFaceColor = [0 0.15 0.95];
        fdBxChart(i1).BoxWidth = bxWidth;
        fdBxChart(i1).MarkerStyle = '.';
        fdBxChart(i1).MarkerColor = [0 0.15 0.95];
        ylabel(pv{i1},'Interpreter','latex') 
        tmpmnm = min(tmpY);
        tmpmax = max(tmpY);
        if mnm(i1)>tmpmnm
            mnm(i1) = tmpmnm;
        end
        if mxm(i1)<tmpmax
            mxm(i1) = tmpmax;
        end
           

    end

    xticks(1:7)
    xlabel('$N$','Interpreter','latex')
 
    % Figure of loss function
    figure(2)
    fdBxChartLoss(ifile) = boxchart(tmp*ifile-deviation,  M(:,end)+1e-7); hold on
    fdBxChartLoss(ifile).BoxFaceColor = [0. 0.15 0.95];
    fdBxChartLoss(ifile).BoxWidth = bxWidth;
    fdBxChartLoss(ifile).MarkerStyle = '.';
    fdBxChartLoss(ifile).MarkerColor = [0. 0.15 0.95];
    xticks(1:7)
    xlabel('$N$','Interpreter','latex')
    ylabel('Loss', 'Interpreter','latex')
    box on

end

legend([fdBxChartLoss(i1),inBxChartLoss(i1)],{'Noninvasive','Invasive'})
set(gca,'YScale','log')

% cosmetic figure configuration
figure(1)
for i1 = 1:5
    subplot(3,2,i1)
    plot([0.5 6.5], [0 0]*prm(i1),':' ), hold on
    lm = max(abs([mnm(i1) mxm(i1)]))
    %ylim([-lm lm]*1.2)
    if i1<3
        set(gca,'XTickLabel',[])
    else
        set(gca,'XTickLabel',lbl2)
    end
    xlim([0.5 6.5])
    xticks(1:7)
    legend([fdBxChart(i1),inBxChartPar(i1)],{'Noninvasive','Invasive'})
    box on
end


 for i1 = 1:5
     subplot(3,2,i1)
     set(gca,'Position',ax1(i1,:))    
 end


set(gcf, 'PaperSize', [17 9]);
set(gcf,'Units','centimeters')
set(gcf,'Position',[100,1000,17,9]);

figure(2)
set(gca,'XTickLabel',lbl2)
xlim([0.5 6.5])
xticks(1:7)
set(gcf, 'PaperSize', [10 6]);
set(gcf,'Units','centimeters')
set(gcf,'Position',[100,1000,10,6]);



