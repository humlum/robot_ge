sp.exper=1;
geq = load(char(strcat('../output/geq',filename,exper.title(sp.exper),'.mat')));
figname = char(strcat('geq',filename,exper.title(sp.exper),'_'));

%% Figure 4: Share of Robot Adopters in Manufacturing
fig = figure('Name',char(strcat(figname,'diffusionCurve')));
hold on
plot(env.ySim(1:length(init.adoptionCount)),init.adoptionCount,'o','LineWidth',2.3,'color',figs.color(1,:))
plot(env.ySim,squeeze(sum(geq.frm.density(1,:,2,:),2)),'LineWidth',2.3,'color',figs.color(1,:))
plot(env.ySim,squeeze(sum(geq.frm.density(2,:,2,:),2)),'LineWidth',2.3,'color',figs.color(2,:))
leg = legend({'Data',exper.legend{1,1},exper.legend{1,2}},'location', 'northwest');
leg.FontSize = 15.5;
ylim([0 1]);
set(gca,'FontSize',15)
legend boxoff;
grid on;
hold('off')
output.printfig(fig);

%% Figure 5: Real Wage Effects of Industrial Robots
wagesReal = geq.wages./permute(repmat(geq.cpi,1,1,env.nOcc+1,env.nSectors),[1 3 4 2]);

% Average real wages
wagesRealAvg = squeeze(sum(sum(wagesReal(:,1:end-1,:,:).*permute(repmat(squeeze(sum(sum(sum(geq.wrk.density(1,:,:,:,:,:,figs.tEval),4),3),2)),1,1,2),[3 1 2]),3),2)); % Average real wage effect (fixed shares)
squeeze(log(wagesRealAvg(1,figs.tEval)./wagesRealAvg(2,figs.tEval)))*100

% Production worker real wages
yEval=2019;
occRef=2;secRef=1;
squeeze(log(wagesReal(1,occRef,secRef,yEval-1990+1)./wagesReal(2,occRef,secRef,yEval-1990+1)))*100

% Tech worker real wages
yEval=2019;
occRef=3;secRef=1;
squeeze(log(wagesReal(1,occRef,secRef,yEval-1990+1)./wagesReal(2,occRef,secRef,yEval-1990+1)))*100

% Other occupations
squeeze(log(wagesReal(1,:,:,yEval-1990+1)./wagesReal(2,:,:,yEval-1990+1)))*100

fig = figure('Name',char(strcat(figname,'realWages')));
hold on
i=0;
for s=1:2 
    for o=[3 2 1]
        i=i+1;
            plot(env.ySim,squeeze(log(wagesReal(1,o,s,:)./wagesReal(2,o,s,:)))*100,'Color',figs.colorOrder(i,:),'LineWidth',2.5)
    end
end
grid on;
ylabel('Log points (percent)','FontSize',15)
y1=get(gca,'ylim');
set(gca,'FontSize',15)
min1 = y1(1) - 0.05*(y1(2)-y1(1));
max1 = y1(2) + 0.05*(y1(2)-y1(1));
plot([2019 2019],[min1 max1],'--k');
ylim([min1 max1]);
legend(figs.legendOrder,'location', 'southwest','FontSize', 11.5);
legend boxoff;
hold('off')
output.printfig(fig);



%% Figure 6: Decomposition of the Production Worker Wage Effect
supplyLabor = nan(2,env.nOcc,env.nSectors,env.nYears);
demand = nan(2,env.nSectors,env.nYears);
factorShareManuf = nan(2,env.nOcc+1,env.nYears);
wbAgg = nan(2,env.nOcc,env.nSectors,env.nYears);

% Calculate labor supply + Demand
for costNo=1:2
    supplyLabor(costNo,:,:,:) = permute(repmat(init.wrk.mass,env.nOcc,1,env.nSectors),[1 3 2]).*squeeze(sum(sum(sum(squeeze(geq.wrk.density(costNo,:,:,:,:,:,:,:)).*par.wrk.hcap,3),2),1));
    demand(costNo,:,:) = func.demand(env,par,squeeze(geq.frm.density(costNo,:,:,:)),squeeze(supplyLabor(costNo,:,:,:)),squeeze(geq.wages(costNo,:,:,:)));
    factorShareManuf(costNo,:,:) = func.AggFactorShare(env,par,squeeze(geq.frm.density(costNo,:,:,:)),squeeze(geq.wages(costNo,:,1,:)));
    wbAgg(costNo,:,2,:) = squeeze(par.alpha(1:end-1,2,:)).*squeeze(demand(costNo,2,:)).';
    wbAgg(costNo,:,1,:) = (squeeze(factorShareManuf(costNo,1:end-1,:)).*squeeze(demand(costNo,1,:)).')./par.frm.markup;
end

% Calculate labor supply + demand + CPI effects
occRef = 1;
secRef = 2;
dl = struct(); % Effect of robots (logarithmic differences)
dl.wages = squeeze(log(geq.wages(1,1:end-1,:,:))-log(geq.wages(2,1:end-1,:,:)));
dl.supplyLabor = squeeze(log(supplyLabor(1,:,:,:)./supplyLabor(1,occRef,secRef,:))-log(supplyLabor(2,:,:,:)./supplyLabor(2,occRef,secRef,:))); % laborDemand Effects are relative to other workers in services
dl.laborDemand = squeeze(log(wbAgg(1,:,:,:)./wbAgg(1,occRef,secRef,:))-log(wbAgg(2,:,:,:)./wbAgg(2,occRef,secRef,:))); % laborDemand Effects are relative to other workers
dl.cpi = squeeze(log(geq.cpi(1,:))-log(geq.cpi(2,:)));
dl.realWages = dl.wages - permute(repmat(dl.cpi,env.nOcc,1,env.nSectors),[1 3 2]);

% Wage decomposition by occ 
for occ=2
    for sec=1
        fig = figure('Name',char(strcat(figname,sprintf('wageDecomposition_occ%d_sec%d',occ,sec))));
        hold on
        plot(env.ySim,squeeze(dl.realWages(occ,sec,:))*100,'LineWidth',2.5)
        plot(env.ySim,squeeze(dl.laborDemand(occ,sec,:))*100,'--','LineWidth',2.3)
        plot(env.ySim,squeeze(-dl.supplyLabor(occ,sec,:))*100,'--','LineWidth',2.3)
        plot(env.ySim,squeeze(-dl.cpi)*100,'--','LineWidth',2.3)
        ylabel('Log Points  (percent)')
        legend({'Real Wages','Labor Demand','Labor Supply','Consumer Prices'},'location', 'northwest','FontSize', 14.5)
        hold off
        ylabel('Log points (percent)','FontSize',15)
        set(gca,'FontSize',15)
        ylim([-50 50]);
        yticks(-40:20:40)
        legend boxoff;
        grid on;
        output.printfig(fig);
    end
end


%% Figure 7: Welfare Effects for Manufacturing Production Workers in 2019
% Welfare Effect in 2019: Compensating Variations as in Caliendo, Dvorkin & Parro (2019)
compVar = nan(size(squeeze(geq.wrk.v(1,:,:,:,:,:,1))));
for a=1:env.wrk.nAge
     compVar(:,a,:,:,:) = squeeze(geq.wrk.v(1,:,a,:,:,:,figs.tEval)-geq.wrk.v(2,:,a,:,:,:,figs.tEval))*((1-par.beta)./(1-par.beta^(env.wrk.nAge-a+1))).*100;
end
compVar(isnan(compVar))=0;

% Calculate discounted log life real earnings for workers that stay in occupation for remainder of career
earnings = squeeze(permute(repmat(wagesReal(:,1:end-1,:,:),1,1,1,1,env.wrk.nSkills,env.wrk.nAge,env.wrk.nTen),[1 5 6 7 2 3 4]).*permute(repmat(par.wrk.hcap,1,1,1,1,1,1,2),[7 1 2 3 4 5 6]));
wOcc = zeros(2,env.wrk.nSkills,env.wrk.nAge,env.wrk.nTen,env.nOcc,env.nSectors);
for costNo=1:2
    for a=1:env.wrk.nAge
        for t=1:min(a,env.wrk.nTen)
            for z=1:env.wrk.nAge-a+1
                wOcc(costNo,:,a,t,:,:) = wOcc(costNo,:,a,t,:,:) + par.beta^(z-1).*func.util(earnings(costNo,:,min(a+z-1,env.wrk.nAge),min(t+z-1,env.wrk.nTen),:,:,min(figs.tEval+z-1,env.nYears)));      
            end
        end
        wOcc(costNo,:,a,:,:,:) = wOcc(costNo,:,a,:,:,:)./((1-par.beta.^(env.wrk.nAge-a+1))./(1-par.beta));
    end  
end
dl.wOcc = squeeze(wOcc(1,:,:,:,:,:) - wOcc(2,:,:,:,:,:))*100;

% Manufacturing Production workers (low-skilled): Welfare Effects Decomposition By Age 
skill = 1;
occ = 2;
sec = 1;

dl.wOccAge = squeeze(sum(squeeze(dl.wOcc(1,:,:,:,:)).*squeeze(geq.wrk.density(costNo,1,:,:,:,:,figs.tEval)),2)./sum(squeeze(geq.wrk.density(costNo,1,:,:,:,:,figs.tEval)),2));
dl.valueAge = squeeze(sum(squeeze(compVar(1,:,:,:,:)).*squeeze(geq.wrk.density(costNo,1,:,:,:,:,figs.tEval)),2)./sum(squeeze(geq.wrk.density(costNo,1,:,:,:,:,figs.tEval)),2));
dl.optionAge = dl.valueAge-dl.wOccAge;

aGrid = linspace(25,65,(env.wrk.nAge+1)/5+1);
xGrid = linspace(1,40,(env.wrk.nAge+1)/5+1); %-0.5;

fig = figure('Name',char(strcat(figname,sprintf('welfare_optionvalue_occ%d_sec%d_y%d',occ,sec,env.ySim(figs.tEval)))));
hold on
bar(squeeze(dl.wOccAge(:,occ,sec)),'FaceColor',figs.color(2,:),'FaceAlpha',0.5,'EdgeColor','none');
bar(squeeze(dl.optionAge(:,occ,sec)),'FaceColor',figs.color(1,:),'FaceAlpha',0.5,'EdgeColor','none');
plot(squeeze(dl.valueAge(:,occ,sec)),'x','color',[0 0 0],'MarkerSize',11,'LineWidth',2.3);
set(gca,'XTick',xGrid);
xlim([0 40.5]);
xlabel('Age','FontSize',15)
ylabel('Percent of Remaining Lifetime Earnings','FontSize',15)
legend({'Production Wages','Option Value','Welfare'},'FontSize',14)
ylabel('Log Points (percent)')
grid on;
set(gca,'XTickLabel',aGrid);
hold off
ylim([-12 12])
set(gca,'FontSize',15)
legend boxoff;
output.printfig(fig);


%% Figure 8: The Effect of Industrial Robots on Employment Shares
for sec=1 
    for occ=2:3 
        fig = figure('Name',char(strcat(figname,sprintf('supplyLabor_sec%d_occ%d',sec,occ))));
        hold on
        plot(env.ySim(1:length(init.adoptionCount)),squeeze(init.supplyLabor(occ,sec,1:length(init.adoptionCount)))*100,'O','LineWidth',2.5)
        plot(env.ySim,squeeze(supplyLabor(1,occ,sec,:))*100,'LineWidth',2.5,'Color',figs.color(1,:))
        plot(env.ySim,squeeze(supplyLabor(2,occ,sec,:))*100,'LineWidth',2.5,'Color',figs.color(2,:))
        if occ==2
            legend({'Data',exper.legend{1,1},exper.legend{1,2}},'location', 'southwest','FontSize',17); 
        else
            legend({'Data',exper.legend{1,1},exper.legend{1,2}},'location', 'southeast','FontSize',17); 
        end
        ylabel('Percent','FontSize',16);
        set(gca,'FontSize',16)
        legend boxoff;
        ylim([1.3 4.7])
        grid on;
        hold('off')
        output.printfig(fig);
    end
end

% Employment share changes accounted for by robotization
disp('Change in Employment Shares: Data')
squeeze(supplyLabor(1,:,:,figs.tEval)-supplyLabor(1,:,:,1))*100
disp('Effect of Robots (percentage points)')
squeeze(supplyLabor(1,:,:,figs.tEval)-supplyLabor(2,:,:,figs.tEval))*100
disp('Share of data accounted for by robots')
squeeze(supplyLabor(1,:,:,figs.tEval)-supplyLabor(2,:,:,figs.tEval))./squeeze(supplyLabor(1,:,:,figs.tEval)-supplyLabor(1,:,:,1))*100