classdef output % Auxiliary functions for analyzing GE output
    methods(Static)
        function [] = printfig(figin)
            fig = figure(figin);
            fig.PaperUnits = 'centimeters';
            fig.PaperPositionMode = 'manual';
            fig.PaperPosition = [0 0 16 12];
            fig.PaperSize = [16 12];
            filename = ['figs\' get(fig,'name') ''];
            print('-dpdf',['' filename '.pdf']);
            
        end
        
        
        % Calculate endogenous GE variables
        function geq = geVar(env,par,init,sol,geq)
            
            % Arrays for endogenous variables
            % General economy
            geq.supplySkills = nan(2,env.nOcc,env.nSectors,env.nYears);
            geq.cpi = nan(2,env.nYears);
            
            % Workers
            geq.wrk.v = nan(2,env.wrk.nSkills,env.wrk.nAge,env.wrk.nTen,env.nOcc,env.nSectors,env.nYears);
            
            % Firms
            geq.frm.sales = nan(2,par.frm.zPoints,2,env.nYears);
            geq.frm.profit = nan(2,par.frm.zPoints,2,env.nYears);
            geq.frm.demand = nan(2,env.nSectors,env.nYears);
            geq.frm.v = nan(2,par.frm.zPoints,2,env.nYears);
            geq.frm.adopt = nan(2,par.frm.zPoints,env.nYears);
            
            % Compute endogenous variables
            for costNo=1:2
                % General Economy
                geq.supplySkills(costNo,:,:,:) = permute(repmat(init.wrk.mass,env.nOcc,1,env.nSectors),[1 3 2]).*squeeze(sum(sum(sum(squeeze(geq.wrk.density(costNo,:,:,:,:,:,:)).*par.wrk.hcap,3),2),1));
                geq.cpi(costNo,:) = func.cpi(env,par,squeeze(geq.wages(costNo,:,:,:)),squeeze(geq.frm.density(costNo,:,:,:))); % cpi
                
                % Workers
                [geq.wrk.v(costNo,:,:,:,:,:,:), ~, ~] = solve.wrk(env,par,squeeze(geq.wages(costNo,:,:,:))./permute(repmat(geq.cpi(costNo,:),env.nOcc+1,1,env.nSectors),[1 3 2]),1); % geq.wrk.v
                
                % Firms
                [geq.frm.sales(costNo,:,:,:), geq.frm.profit(costNo,:,:,:)] = func.frmFlow(env,par,squeeze(geq.wages(costNo,:,:,:)),squeeze(geq.frm.density(costNo,:,:,:)),squeeze(geq.supplySkills(costNo,:,:,:))); % msm.sales
                geq.frm.demand(costNo,:,:) = func.demand(env,par,squeeze(geq.frm.density(costNo,:,:,:)),squeeze(geq.supplySkills(costNo,:,:,:)),squeeze(geq.wages(costNo,:,:,:))); % geq.frm.demand
                geq.frm.v(costNo,:,:,:) = solve.frm(env,par,sol,squeeze(geq.frm.profit(costNo,:,:,:))./permute(repmat(geq.cpi(costNo,:),2,1,par.frm.zPoints), [3 1 2]),geq.cRobot(costNo,:)./geq.cpi(costNo,:),1);
                [~, geq.frm.adopt(costNo,:,:)] = simulate.frm(env,par,squeeze(geq.frm.density(costNo,:,:,:)),geq.cRobot(costNo,:)./geq.cpi(costNo,:),squeeze(geq.frm.v(costNo,:,:,:)),1);
                
            end
        end
        
        function [cvCurrent, cvCurrentProd, cvFuture] = cvTax(env,par,init,geq,get,tEval)
            wagesReal = geq.wages./permute(repmat(geq.cpi,1,1,env.nOcc+1,env.nSectors),[1 3 4 2]);
            earnings = permute(repmat(wagesReal(1,1:end-1,:,:),env.wrk.nSkills,1,1,1,env.wrk.nAge,env.wrk.nTen),[1 5 6 2 3 4]).*par.wrk.hcap; % Compensating variation earnings depends on current income            
            compVar = zeros(size(squeeze(get.wrk.v(1,:,:,:,:,:,1))));
            cvIncome = zeros(size(compVar));
            for a=1:env.wrk.nAge
                compVar(:,a,:,:,:) = exp(squeeze(get.wrk.v(2,:,a,:,:,:,tEval)-get.wrk.v(1,:,a,:,:,:,tEval))*((1-par.beta)./(1-par.beta^(env.wrk.nAge-a+1))))-1;
                cvIncome(:,a,:,:,:) = squeeze(compVar(:,a,:,:,:)).*squeeze(earnings(:,a,:,:,:,tEval)).*((1-par.beta^(env.wrk.nAge-a+1))./(1-par.beta));         
            end
            cvIncome(isnan(cvIncome))=0;
            cvCurr = init.wrk.mass(tEval).*sum(sum(sum(sum(sum(cvIncome.*squeeze(get.wrk.density(2,:,:,:,:,:,tEval)))))));
            cvCurrent = init.wrk.mass(tEval).*sum(cvCurr(:));
            cvCurrentProd = init.wrk.mass(tEval).*sum(sum(sum(sum(sum(cvIncome(:,:,:,2,1).*squeeze(get.wrk.density(2,:,:,:,2,1,tEval)))))));
            
            compVar = nan(env.wrk.nSkills,env.nOcc,env.nSectors);
            for t=tEval+1:env.nYears
                compVar = exp(squeeze(get.wrk.v(2,:,1,1,:,:,t)-get.wrk.v(1,:,1,1,:,:,t))*((1-par.beta)./(1-par.beta^(env.wrk.nAge))))-1;
                cvFuture(t) = par.beta^(t-tEval).*init.wrk.mass(t).*sum(sum(sum(sum(squeeze(get.wrk.density(2,:,1,1,:,:,t)).*compVar.*squeeze(earnings(:,1,1,:,:,t)).*((1-par.beta^(env.wrk.nAge))./(1-par.beta)))))); % .*((((1+gRate)^(env.wrk.nAge)-1)/gRate)/(env.wrk.nAge));
            end
            cvFuture(:,end) = cvFuture(:,end)./(1-par.beta);
            cvFuture = sum(cvFuture,2);
        end
        
        function [taxPDV, taxPDVInelastic] = taxPDV(env,par,get,tEval,permanent,tFullAdopt)
            taxPDV = 0;
            taxPDVInelastic = 0;
            for t=tEval:env.nYears-1
                taxFlow =            squeeze(get.cRobot(2,t)-get.cRobot(1,t))*sum(squeeze(get.frm.adopt(2,:,t)).*squeeze(get.frm.density(2,:,1,t)),2)./squeeze(get.cpi(2,t));
                taxPDV = taxPDV+par.beta^(t-tEval)*taxFlow;
                taxFlowInelastic =   squeeze(get.cRobot(2,t)-get.cRobot(1,t))*sum(squeeze(get.frm.adopt(1,:,t)).*squeeze(get.frm.density(1,:,1,t)),2)./squeeze(get.cpi(2,t));
                taxPDVInelastic = taxPDVInelastic+par.beta^(t-tEval)*taxFlowInelastic;
                if permanent==1 && t==env.nYears-1
                    taxPDV =            taxPDV          + par.beta^(env.nYears-tEval)*squeeze(get.cRobot(2,end)-get.cRobot(1,end))*((1-sum(get.frm.density(2,:,2,end),2))/tFullAdopt)*((1-par.beta^(tFullAdopt+1))/(1-par.beta))./squeeze(get.cpi(2,t));
                    taxPDVInelastic =   taxPDVInelastic + par.beta^(env.nYears-tEval)*squeeze(get.cRobot(2,end)-get.cRobot(1,end))*((1-sum(get.frm.density(1,:,2,end),2))/tFullAdopt)*((1-par.beta^(tFullAdopt+1))/(1-par.beta))./squeeze(get.cpi(2,t));
                end
            end
        end
        
        function profitPDV = profitPE(env,sol,par,init,get,tEval)
            frmV = nan(size(get.frm.v));
            for costNo=1:2
                cost = permute(repmat(par.costAdj(1,:),par.frm.zPoints,1,2),[1 3 2]).*func.cost(env,par,squeeze(get.wages(costNo,:,1,:)));
                price = func.priceIndex(par,cost,init.frm.density); % This is the restriction that shuts off product market stealing effects: init.frm.density
                
                supplySkills = permute(repmat(init.wrk.mass,env.nOcc,1,env.nSectors),[1 3 2]).*squeeze(sum(sum(sum(squeeze(get.wrk.density(costNo,:,:,:,:,:,:)).*par.wrk.hcap,3),2),1));
                demand = func.demand(env,par,init.frm.density,supplySkills,squeeze(get.wages(costNo,:,:,:)));
                profit = func.profitGrid(par,cost,price,demand(1,:));
                
                frmV(costNo,:,:,:) = solve.frm(env,par,sol,profit./permute(repmat(get.cpi(costNo,:),2,1,par.frm.zPoints), [3 1 2]),get.cRobot(costNo,:)./get.cpi(costNo,:),1);
            end
            profitPDV = sum(sum(squeeze((frmV(2,:,:,tEval)-frmV(1,:,:,tEval)).*get.frm.density(1,:,:,tEval)),2),1);
        end
        
        function incidence = incidence(get,gdp)
            incidence = nan(8,1);
            incidence(1,:) = round((get.cvCurrent+get.cvFuture)/gdp*100,2).';
            incidence(2,:) = round(get.cvCurrent/gdp*100,2).';
            incidence(3,:) = round(get.cvProd./gdp*100,2).';
            incidence(4,:) = round(get.cvFuture/gdp*100,2).';
            incidence(5,:) = round(get.taxPDV/gdp*100,2).';
            incidence(6,:) = round(get.taxPDVInelastic/gdp*100,2).';
            incidence(7,:) = round((get.taxPDV-get.taxPDVInelastic)/gdp*100,2).';
            incidence(8,:) = round(get.profitPDV/gdp*100,2).';
%            incidence = table(incidence);
%            
%            incidenceCol = [get.legend; num2cell(incidence)];
        end
        
        function incidenceTable = incidenceTable(get1,get2,gdp)
            get1.incidence = output.incidence(get1,gdp);
            get2.incidence = output.incidence(get2,gdp);
            incidenceLabel = cell(8,1);            
            incidenceLabel{1,1} = 'Workers';
            incidenceLabel{2,1} = '\hskip0.6cm Workers in 2019';
            incidenceLabel{3,1} = '\hskip0.9cm -- Manufacturing Production';
            incidenceLabel{4,1} = '\hskip0.6cm Future Workers';
            incidenceLabel{5,1} = 'Tax Revenues';
            incidenceLabel{6,1} = '\hskip0.6cm Mechanical Effect';
            incidenceLabel{7,1} = '\hskip0.6cm Behavioral Effect';
            incidenceLabel{8,1} = 'Profits (excl. predatory externalities)';
            incidenceTable = table(get1.incidence,get2.incidence,'RowNames',incidenceLabel,'VariableNames',{char(get1.legend) char(get2.legend)}); % 
 %           incidenceTable.Properties.VariableNames = {char(get1.legend) char(get2.legend)};
 %           horzcat(incidenceLabel,num2cell(get1.incidence),num2cell(get2.incidence)));            
 %           incidenceTable = table(horzcat(incidenceLabel,num2cell(get1.incidence),num2cell(get2.incidence)));
 %            
        end
        
        
    end
end