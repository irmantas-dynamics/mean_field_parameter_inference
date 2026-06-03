close all
clear all



thr = 1e-3;

cmap = makeColorMap([255 255 255]./255,[204 229 148]./255,[10 10 170]./255,128);


figure
load('arnold_QIF_AD.mat','TperM','ampM','dist');
subplot(1,2,2)
imagesc(TperM,ampM,dist<thr), hold on
plot(80,-4,'r+')
colormap(cmap)
set(gca,'YDir','normal')
ylabel('$K$','Interpreter','latex')
xlabel('$T_\mathrm{ext}$, ms','Interpreter','latex')


load('lyap_QIF_AD.mat','KM','maxL');
subplot(1,2,1)
plot(KM,maxL), hold on
plot(KM,KM*0,':')
plot(5,-0.0173,'r+')
ylim([-0.035 0.03])
xlim([-inf 10])
ylabel('$\lambda_\mathrm{max}$, ms$^{-1}$','Interpreter','latex')
xlabel('$K$','Interpreter','latex')


set(gcf, 'PaperSize', [16 4]);
set(gcf,'Units','centimeters')
set(gcf,'Position',[2,2,16,4])  
