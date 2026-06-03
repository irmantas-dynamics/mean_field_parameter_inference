close all
clear all




figure
load("Loss.mat")
subplot(121)
s = imagesc(AMas,rMas,LossM);
xlabel('$A_0$','Interpreter','latex')
ylabel('$R_0$','Interpreter','latex')
set(gca,'ydir','normal')
title("Without synchronization")

clim([0 3.5])

load("Loss_sync.mat")
subplot(122)

s = imagesc(AMas,rMas,LossM);
xlabel('$A_0$','Interpreter','latex')
ylabel('$R_0$','Interpreter','latex')
title("With synchronization")
set(gca,'ydir','normal')
clim([0 3.5])
clb = colorbar;
clb.Label.String = 'Loss value';
hold on

colormap jet


set(gcf, 'PaperSize', [17 9]);
set(gcf,'Units','centimeters')
set(gcf,'Position',[100,1000,17,9]);


