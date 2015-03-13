% =============================
%
% Plotting function for DG
% By Devin Light
%
% =============================

clc;
clear all;
whichMeth = 'variableCoeff';
outDir = ['_output/' whichMeth '/'];
figureDir = '_figures/';

%%
close all;
FS = 'FontSize';

meqn        = 2;

whichRun    = 1;
numFrames   = 11;
plotAxis    = [-5 5 -.2 2.2];

plotCoeffs  = 1;
maxDegree   = 10;
contours    = 1:maxDegree;
printDofs   = 1;

printYlabel = 1;
saveFig     = 0;

xwidth = 400; ywidth = 400;

for n=1:numFrames
    frameid = n-1;
    [out aux] = readData(outDir,whichRun,frameid);
%    fig = figure();
%    set(gcf, 'PaperUnits', 'points');
%    set(gcf,'PaperPositionMode','auto','PaperSize',[xwidth ywidth]);
%    set(fig, 'Position', [0 0 xwidth ywidth]);

    x = out.x; t = out.t; data = out.data;
    for m=1:meqn
        q = data(:,m);
        %plot(x,q,'r-');
        title(['q',num2str(m),' at time T=',num2str(t)],FS,16);

        elemBdrys = aux.elemCent';
        dx = elemBdrys(2)-elemBdrys(1);
        elemBdrys = [aux.elemCent(1)-0.5*dx, elemBdrys + 0.5*dx];

        set(gca,'XTick',elemBdrys,'XTickLabel','');
        axis(plotAxis);
        xlabel('x',FS,18); ylabel('y',FS,18);
        set(gca,FS,16);
        
        if(plotCoeffs == 1)
        % Make "heat plot" for coefficient values
            fig = figure();
            set(gcf, 'PaperUnits', 'points');
            set(gcf,'PaperPositionMode','auto','PaperSize',[xwidth ywidth]);
            set(fig, 'Position', [0 0 xwidth ywidth])

            y = linspace(plotAxis(3),plotAxis(4),10);
            phi = zeros(length(y),length(x));
            for j=1:length(aux.elemCent)
                xloc = elemBdrys(j) <= x & x <= elemBdrys(j+1);
                phi(:,xloc) = aux.localPolyDegree(j);
            end

            cMap = colormap(hsv(maxDegree));
            contourf(x,y,phi,contours,'LineWidth',0.5); caxis([0,maxDegree]);
            hold on, plot(x,q,'k-');hold off; axis(plotAxis);
            xlabel('x',FS,18);set(gca,FS,16);
            set(gca,'XTick',elemBdrys,'XTickLabel','');

            foo = phi-phi(1);
            if(sum(foo(:)) == 0)
                % If phi is constant contourf won't plot anything..
                % so change background color to the color for that (constant)
                % value
                set(gcf,'Color','white');
                set(gca,'Color',cMap(phi(1),:));
                set(gcf, 'InvertHardCopy', 'off');
            end
            
            if(printYlabel == 1)
                ylabel('y',FS,18); 
            else
                ylabel('');set(gca,'YTickLabel','');
            end

            if(printDofs == 1)
                % Print number of dofs in top right corner
                xloc = 0.7*(plotAxis(2)-plotAxis(1));
                yloc = 0.8*(plotAxis(4)-plotAxis(3));
                dofCount = sum(aux.localPolyDegree);
                text(xloc,yloc,['DOF = ' sprintf('%4i',dofCount)],FS,16);
            end

            title(['q',num2str(m),' at time T=',num2str(t)],FS,16);
            name = [figureDir whichMeth '_q' num2str(m) 'coeffContour_run' num2str(whichRun) '_frame' num2str(n) '.pdf'];

            if(saveFig == 1)
                print(fig,name,'-dpdf');
            end
            
        end
    end
end

if(plotCoeffs == 1)
    % Print colorbar figure for coeff heat plot
    fig = figure(); axis off;
    set(gcf, 'PaperUnits', 'points');
    set(gcf,'PaperPositionMode','auto','PaperSize',[xwidth ywidth]);
    set(fig, 'Position', [0 0 xwidth ywidth])
    colormap(cMap); colorbar('location','eastoutside',FS,16); 
    caxis([0,maxDegree]);
    name = [figureDir 'coeffContour_run' num2str(whichRun) '_cb.pdf'];
    print(fig,name,'-dpdf');
end

pause(0.2); %close all;

%% Make efficiency plot
MS = 'MarkerSize';

maxDegree = 10;
ne = 16*[1 2 4 8 16];
unAdaptDOF = (maxDegree+1)*ne;
adaptDOF = 0.*unAdaptDOF;

xwidth = 400; ywidth = 400;
fig = figure(); axis off;
set(gcf, 'PaperUnits', 'points');
set(gcf,'PaperPositionMode','auto','PaperSize',[xwidth ywidth]);
set(fig, 'Position', [0 0 xwidth ywidth])
    
outDir = '_output/padg/';
for n=1:length(ne)
    [out aux] = readData(outDir,n,10);
    adaptDOF(n) = sum(aux.localPolyDegree(:));
end

adaptL2 = [0.6979E-04 0.1484E-04 0.2799E-05 0.5548E-06 0.1127E-06];
unadaptL2 = [0.7640E-04 0.1658E-04 0.3815E-05 0.9081E-06 0.2439E-06];

appxOrder = polyfit(log10(adaptDOF),log10(adaptL2),1);
tmp = ['   Approximate order = ', num2str(appxOrder(1))];
orderLine = 10.^(appxOrder(2)).*10.^(appxOrder(1)*log10(adaptDOF));

plot(adaptDOF,adaptL2,'k.-',unAdaptDOF,unadaptL2,'r.-','markers',35); 
set(gca,'xscale','log','yscale','log');
xlabel('DOF',FS,16); ylabel('L2 Error',FS,16);
axis([56.864170679821,5686.417067982086,0.000000024696649422,0.000246966494218378]);
leg=legend('p-Adaptive','non-adaptive');
set(leg,FS,16); set(gca,FS,14);

name = 'efficency.pdf';
print(fig,name,'-dpdf');





