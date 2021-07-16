geq = load(char(strcat('../output/geq',filename,exper.title(1),'.mat')));
geqT = load(char(strcat('../output/geq',filename,exper.title(2),'.mat')));
geqP = load(char(strcat('../output/geq',filename,exper.title(3),'.mat')));
figname = char(strcat('geq',filename,exper.title(2),'_',exper.title(3),'_'));
geqT.legend = exper.legend(2,2);
geqP.legend = exper.legend(3,2);

%% Figure 9: Robot Tax Counterfactuals
tEnd=1;
densityAdopt = (squeeze(squeeze(geq.frm.density(1,:,1,1:tEnd))).*squeeze(geq.frm.adopt(1,:,1:tEnd)))./sum(squeeze(squeeze(geq.frm.density(1,:,1,1:tEnd))).*geq.frm.adopt(1,:,1:tEnd),2);
salesAdopt = sum(squeeze(squeeze(geq.frm.sales(1,:,1,1:tEnd)).*densityAdopt(1,:)),2); 

% Figure 9(a): Robot Adoption Costs
fig = figure('Name',char(strcat(figname,'cRobot')));
hold on
plot(env.ySim,squeeze(geqT.cRobot(1,:))./salesAdopt,'-','LineWidth',2.5,'Color',figs.color(1,:))
plot(env.ySim,squeeze(geqT.cRobot(2,:))./salesAdopt,'LineWidth',2.5,'Color',figs.color(2,:))
plot(env.ySim,squeeze(geqP.cRobot(2,:))./salesAdopt,'LineWidth',2.5,'Color',figs.color(3,:))
y1=get(gca,'ylim');
max1 = y1(2) + 0.05*(y1(2)-y1(1));
set(gca,'FontSize',16)
ylabel('Units of adopter sales','FontSize',16)
ylim([0 max1]);
plot([env.ySim(exper.tSurprise(2)) env.ySim(exper.tSurprise(2))],[0 max1],'--k');
plot([env.ySim(exper.tSurprise(2)+10) env.ySim(exper.tSurprise(2)+10)],[0 max1],'--k');
legend({exper.legend{2,1},exper.legend{2,2},exper.legend{3,2}},'location', 'southwest','FontSize',17);
grid on;
legend boxoff;
hold('off')
figures.printfig(fig);

%% Figure 9(b): Robot Diffusion Curve
fig = figure('Name',char(strcat(figname,'diffusionCurve')));
hold on
plot(env.ySim(1:length(init.adoptionCount)),init.adoptionCount,'o','LineWidth',2.5)
plot(env.ySim,squeeze(sum(geqT.frm.density(1,:,2,:),2)),'LineWidth',2.5,'Color',figs.color(1,:))
plot(env.ySim,squeeze(sum(geqT.frm.density(2,:,2,:),2)),'LineWidth',2.5,'Color',figs.color(2,:))
plot(env.ySim,squeeze(sum(geqP.frm.density(2,:,2,:),2)),'LineWidth',2.5,'Color',figs.color(3,:))
ylabel('Share of robot adopters in manufacturing','FontSize',16);
ylim([0 1]);
plot([env.ySim(exper.tSurprise(2)) env.ySim(exper.tSurprise(2))],[0 1],'--k');
plot([env.ySim(exper.tSurprise(2)+10) env.ySim(exper.tSurprise(2)+10)],[0 1],'--k');
set(gca,'FontSize',16)
legend({'Data',exper.legend{2,1},exper.legend{2,2},exper.legend{3,2}},'location', 'northwest','FontSize',17); %northwest southeast
legend boxoff;
grid on;
hold('off')
figures.printfig(fig);

%% Table 3: Robot Tax Incidence
% GDP in 2019
wrkIncome = init.wrk.mass(1,figs.tEval).*sum(sum(squeeze(sum(sum(sum(squeeze(geq.wrk.density(1,:,:,:,:,:,figs.tEval)).*par.wrk.hcap(:,:,:,:,:,figs.tEval),3),2),1)).*squeeze(geq.wages(1,1:end-1,:,figs.tEval)),2),1)./geq.cpi(1,figs.tEval);
frmProfit = (geq.frm.demand(1,1,figs.tEval)/par.frm.epsilon)./geq.cpi(1,figs.tEval);
gdp = wrkIncome+frmProfit;
wrkIncome/gdp

% Workers 
[geqT.cvCurrent, geqT.cvProd, geqT.cvFuture] = output.cvTax(env,par,init,geq,geqT,figs.tEval);
[geqP.cvCurrent, geqP.cvProd, geqP.cvFuture] = output.cvTax(env,par,init,geq,geqP,figs.tEval);

% Tax revenues 
[geqT.taxPDV, geqT.taxPDVInelastic] = output.taxPDV(env,par,geqT,figs.tEval,0,0);
[geqP.taxPDV, geqP.taxPDVInelastic] = output.taxPDV(env,par,geqP,figs.tEval,1,20);

% Profits (without Product Market Stealing): Take CPI into account still?
geqT.profitPDV = output.profitPE(env,sol,par,init,geqT,figs.tEval);
geqP.profitPDV = output.profitPE(env,sol,par,init,geqP,exper.tSurprise(3));

% Table
incidenceTab = output.incidenceTable(geqT,geqP,gdp);
table2latex(incidenceTab,char(strcat('../output/',figname,sprintf('incidenceTab'))));