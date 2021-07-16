classdef func
    methods(Static)
        function cost = cost(env,par,wages) % funFirmCES.cost
            costMatrix1 = exp(-(repmat(par.frm.zGrid.',1,2,env.nYears) + permute(repmat([0 1].',1,par.frm.zPoints,env.nYears),[2 1 3]).*permute(repmat(par.frm.gammaProd,par.frm.zPoints,1,2),[1 3 2])));
            costMatrix2 = (sum(exp(permute(repmat(par.frm.phiOcc,1,1,2),[1 3 2]) + permute(repmat([0 1].',1,env.nOcc,env.nYears),[2 1 3]).*repmat(par.frm.gammaOcc,1,2,env.nYears)).*permute(repmat(wages(1:end-1,:),1,1,2),[1 3 2]).^(1-par.frm.sigma),1) + permute(repmat(wages(end,:),1,1,2),[1 3 2]).^(1-par.frm.sigma)).^(1/(1-par.frm.sigma));
            cost = costMatrix1.*repmat(costMatrix2,par.frm.zPoints,1,1);
        end
        
        function priceIndex = priceIndex(par,costGrid,density) % robotModel.priceIndex
            % F is a Nx2 the distribution of idiosyncratic states (R,z)
            % L is aggregate labor supply
            priceIndex = squeeze(sum(dot((par.frm.markup.*costGrid).^(1-par.frm.epsilon),density,2),1).^(1/(1-par.frm.epsilon))).';
        end
        
        
        function sRobot = sRobot(par,cost,density) % robotModel.sRobot
            priceIndex = func.priceIndex(par,cost,density);
            sRobot = squeeze(sum(sum(((par.frm.markup.*cost(:,2,:)).^(1-par.frm.epsilon)).*density(:,2,:),2),1)).'./(priceIndex.^(1-par.frm.epsilon));
        end
        
        function factorShare = factorShare(env,par,wages) % funFirmCES.factorShare
            factorShare = nan(2,env.nOcc+1,env.nYears);
            for r=0:1
                factorShare(r+1,1:end-1,:) = (exp(par.frm.phiOcc+r.*par.frm.gammaOcc).*(wages(1:end-1,:).^(1-par.frm.sigma)))./(sum(exp(par.frm.phiOcc+r.*par.frm.gammaOcc).*(wages(1:end-1,:)).^(1-par.frm.sigma),1)+wages(end,:).^(1-par.frm.sigma));
                factorShare(r+1,end,:) = (wages(end,:).^(1-par.frm.sigma))./(sum(exp(par.frm.phiOcc+r.*par.frm.gammaOcc).*(wages(1:end-1,:)).^(1-par.frm.sigma),1)+wages(end,:).^(1-par.frm.sigma));
            end
        end
        
        function AggFactorShare = AggFactorShare(env,par,density,wages) % funFirmCES.AggFactorShare
            cost = func.cost(env,par,wages);
            factorShare = func.factorShare(env,par,wages);
            sRobot =  func.sRobot(par,cost,density);
            AggFactorShare = sRobot.*squeeze(factorShare(2,:,:)) + (1-sRobot).*squeeze(factorShare(1,:,:));
        end
        
        function demand = demand(env,par,frmDensity,supplyLabor,wages) % segModel.demand
            factorShare = func.AggFactorShare(env,par,frmDensity,squeeze(wages(:,1,:))); %AggFactorShare(env,par,density,wages)
            leontief = nan(env.nSectors,env.nSectors,env.nYears);
            leontiefInv = nan(env.nSectors,env.nSectors,env.nYears);
            for t=1:env.nYears
                leontief(:,:,t) = eye(env.nSectors) - [(par.mu)*factorShare(end,t) (par.mu)*par.alpha(end,2,t); (1-par.mu)*factorShare(end,t) (1-par.mu)*par.alpha(end,2,t)];
                leontiefInv(:,:,t) = inv(leontief(:,:,t));
            end
            income = (par.mu.*par.markup(1) + (1-par.mu).*par.markup(2)).*sum(sum(wages(1:end-1,:,:).*supplyLabor,2),1);
            demand = nan(env.nSectors,env.nYears);
            for t=1:env.nYears
                demand(:,t) = leontiefInv(:,:,t)*[par.mu; 1-par.mu].*income(t);
            end
        end
        
        function demand = demandFrmRes(env,par,frmDensity,supplyLabor,wages,cRobot,adopt) % func.demand but with profits and robot adoption costs (deterministic component) in resource constraint
            factorShare = func.AggFactorShare(env,par,frmDensity,squeeze(wages(:,1,:))); %AggFactorShare(env,par,density,wages)
            leontief = nan(env.nSectors,env.nSectors,env.nYears);
            leontiefInv = nan(env.nSectors,env.nSectors,env.nYears);
            leontiefInvProfit = nan(env.nSectors,env.nSectors,env.nYears);
            costAdopt = cRobot.*sum(squeeze(frmDensity(:,1,:)).*adopt,1);
            for t=1:env.nYears
                leontief(:,:,t) = eye(env.nSectors) - [(par.mu)*factorShare(end,t) (par.mu)*par.alpha(end,2,t); (1-par.mu)*factorShare(end,t) (1-par.mu)*par.alpha(end,2,t)];
                leontiefInv(:,:,t) = inv(leontief(:,:,t));
            end
            income = sum(sum(wages(1:end-1,:,:).*supplyLabor,2),1) + costAdopt;
            demand = nan(env.nSectors,env.nYears);
            for t=1:env.nYears
                leontiefInvProfit(:,t) = inv(eye(env.nSectors)-squeeze(leontiefInv(:,:,t))*[par.mu; 1-par.mu]*[1/par.frm.epsilon 0])*squeeze(leontiefInv(:,:,t)); %#ok<MINV>
                demand(:,t) = leontiefInvProfit(:,:,t)*[par.mu; 1-par.mu].*income(t);
            end
        end
        
        
        function profitGrid = profitGrid(par,costGrid,priceIndex,expend) % robotModel.profitGrid
            % F is a Nx2 the distribution of idiosyncratic states (R,z)
            % L is aggregate labor supply
            prices = permute(repmat(priceIndex,par.frm.zPoints,1,2),[1 3 2]);
            income = permute(repmat(expend,par.frm.zPoints,1,2),[1 3 2]);
            profitGrid = costGrid.^(1-par.frm.epsilon).*prices.^(par.frm.epsilon-1).*income.*(((par.frm.epsilon-1)^(par.frm.epsilon-1))/(par.frm.epsilon^(par.frm.epsilon)));
        end
        
        function u = util(x) % model.util
            u = log(x);
        end
        
        function u = uFrm(x)
            u = x;
%            u = log(x); % NB: Linear flow utility is the reason why larger firms select into robot adoption
        end
        
        function cpi = cpi(env,par,wages,frmDensity)
            frmCost = permute(repmat(par.costAdj(1,:),par.frm.zPoints,1,2),[1 3 2]).*func.cost(env,par,squeeze(wages(:,1,:)));
            frmPrice= func.priceIndex(par,frmCost,frmDensity);
            priceServ = par.costAdj(2,:).*squeeze(prod((wages(:,2,:)./par.alpha(:,2,:)).^(par.alpha(:,2,:)),1)).';
            cpi = ((frmPrice./par.mu).^(par.mu)).*((priceServ./(1-par.mu)).^(1-par.mu));
            cpi = cpi./repmat(cpi(1,1),1,env.nYears);
        end
        
        function [sales, profit] = frmFlow(env,par,wages,frmDensity,supplyLabor)
            demand = func.demand(env,par,frmDensity,supplyLabor,wages);
            frmCost = permute(repmat(par.costAdj(1,:),par.frm.zPoints,1,2),[1 3 2]).*func.cost(env,par,squeeze(wages(:,1,:)));
            frmPrice= func.priceIndex(par,frmCost,frmDensity);
            profit = func.profitGrid(par,frmCost,frmPrice,demand(1,:));
            sales = permute(repmat(demand(1,:),par.frm.zPoints,1,2),[1 3 2]).*((par.frm.markup.*frmCost)./permute(repmat(frmPrice,par.frm.zPoints,1,2),[1 3 2])).^(1-par.frm.epsilon);
        end
        
        function density1 = zLom(par,density0)
           density1 = (density0.'*par.frm.zP).'; 
           % density1 = density0;
        end
        
        
    end
end