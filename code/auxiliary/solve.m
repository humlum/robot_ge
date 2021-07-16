classdef solve
    methods(Static)
        
        %% Solve for firm value functions (solve.frm)
        function v = frm(env,par,sol,profit,cRobot,tStart) %
            v = nan(par.frm.zPoints,2,env.nYears);
            
            % Solve for stationary bellman in R = {0,1}
            w0 = (1/(1-par.beta)).*func.uFrm(squeeze(profit(:,:,end)));
            w1 = w0;
            err=1;
            iter=0;
            while err>sol.tol.frm && iter<sol.iter.frm % sol.iter.frm
                iter = iter+1;
                w1(:,1) = func.uFrm(profit(:,1,end)) + par.frm.zP*par.frm.nu*(log(exp((par.beta.*w0(:,1))/par.frm.nu) + exp((-cRobot(end)+par.beta.*w0(:,2))/par.frm.nu)));
                w1(:,2) = func.uFrm(profit(:,2,end)) + par.beta.*par.frm.zP*(par.frm.theta.*w0(:,1) + (1-par.frm.theta).*w0(:,2));
                err = sum(sum(abs(w1-w0)));
                w0 = 0.5*w0+(1-0.5)*w1;
            end
            v(:,:,end) = w0;
            
            % Backwards induction for t=T-1,...,tStart
            for t=env.nYears-1:-1:tStart
                v(:,1,t) = func.uFrm(profit(:,1,t)) + par.frm.zP*par.frm.nu*(log(exp((-cRobot(t)+par.beta.*v(:,2,t+1))/par.frm.nu) + exp((par.beta.*v(:,1,t+1))/par.frm.nu)));
                v(:,2,t) = func.uFrm(profit(:,2,t)) + par.beta*(par.frm.zP*(par.frm.theta.*v(:,1,t+1) + (1-par.frm.theta).*v(:,2,t+1))); % scale with mean of logit shock (because still receive it (?) even though have adopted)
            end
        end
        
        %% Solve for worker value functions (skill heterogeneity)
        function [v, policy, policyEnter] = wrk(env,par,wages,tStart) %
            eyeArray = zeros(env.wrk.nSkills,env.nOcc,env.nSectors,env.nOcc,env.nSectors);
            eyeArray(:,:,1,:,1) = permute(repmat(eye(env.nOcc),1,1,env.wrk.nSkills),[3 1 2]);
            eyeArray(:,:,2,:,2) = permute(repmat(eye(env.nOcc),1,1,env.wrk.nSkills),[3 1 2]);
            
            tSim = length(wages);
            v = nan(env.wrk.nSkills,env.wrk.nAge,env.wrk.nTen,env.nOcc,env.nSectors,tSim);
            policy = nan(env.wrk.nSkills,env.wrk.nAge,env.wrk.nTen,env.nOcc,env.nSectors,env.nOcc,env.nSectors,tSim);
            policyEnter = nan(env.wrk.nSkills,env.nOcc,env.nSectors,tSim);
            
            % Final period T
            value = nan(env.wrk.nSkills,env.nOcc,env.nSectors,env.nOcc,env.nSectors);
            for a = env.wrk.nAge:-1:1
                for t=1:min(a,env.wrk.nTen)
                    if a == env.wrk.nAge % retiring generation
                        v(:,a,t,:,:,end) = permute(repmat(par.wrk.amenity(:,:,end),1,1,env.wrk.nSkills),[3 1 2 ]) + func.util(permute(repmat(wages(1:end-1,:,end),1,1,env.wrk.nSkills),[3 1 2]).*squeeze(par.wrk.hcap(:,a,t,:,:,end)));
                    elseif a < env.wrk.nAge
                        value = permute(repmat(par.wrk.amenity(:,:,end),1,1,env.nOcc,env.nSectors,env.wrk.nSkills),[5 1 2 3 4]) + func.util(repmat(permute(repmat(wages(1:end-1,:,end),1,1,env.wrk.nSkills),[3 1 2]).*squeeze(par.wrk.hcap(:,a,t,:,:,end)),1,1,1,env.nOcc,env.nSectors)) - squeeze(par.wrk.swCost(:,a,:,:,:,:)) + par.beta.*(eyeArray.*permute(repmat(squeeze(v(:,a+1,min(t+1,env.wrk.nTen),:,:,end)),1,1,1,env.nOcc,env.nSectors),[1 4 5 2 3]) + (1-eyeArray).*permute(repmat(squeeze(v(:,a+1,1,:,:,end)),1,1,1,env.nOcc,env.nSectors),[1 4 5 2 3]));
                        v(:,a,t,:,:,end) = par.wrk.rho.*(double(eulergamma)+log(sum(sum(exp(value./par.wrk.rho),5),4)));
                        policy(:,a,t,:,:,:,:,end) = nan;
                    end
                end
            end
            
            % Backwards induction for t=T-1,T-2,...,tStart
            for y = tSim-1:-1:tStart
                % y
                for a = env.wrk.nAge:-1:1
                    for t=1:min(a,env.wrk.nTen)
                        if a == env.wrk.nAge % retiring generation
                            v(:,a,t,:,:,y) = permute(repmat(par.wrk.amenity(:,:,y),1,1,env.wrk.nSkills),[3 1 2 ]) + func.util(permute(repmat(wages(1:end-1,:,y),1,1,env.wrk.nSkills),[3 1 2]).*squeeze(par.wrk.hcap(:,a,t,:,:,y)));
                        elseif a<env.wrk.nAge
                            value = permute(repmat(par.wrk.amenity(:,:,y),1,1,env.nOcc,env.nSectors,env.wrk.nSkills),[5 1 2 3 4]) + func.util(repmat(permute(repmat(wages(1:end-1,:,y),1,1,env.wrk.nSkills),[3 1 2]).*squeeze(par.wrk.hcap(:,a,t,:,:,y)),1,1,1,env.nOcc,env.nSectors)) - squeeze(par.wrk.swCost(:,a,:,:,:,:)) + par.beta.*(eyeArray.*permute(repmat(squeeze(v(:,a+1,min(t+1,env.wrk.nTen),:,:,y+1)),1,1,1,env.nOcc,env.nSectors),[1 4 5 2 3]) + (1-eyeArray).*permute(repmat(squeeze(v(:,a+1,1,:,:,y+1)),1,1,1,env.nOcc,env.nSectors),[1 4 5 2 3]));
                            v(:,a,t,:,:,y) = par.wrk.rho.*(double(eulergamma)+log(sum(sum(exp(value./par.wrk.rho),5),4)));
                            policy(:,a,t,:,:,:,:,y) = exp(value./par.wrk.rho)./sum(sum(exp(value./par.wrk.rho),5),4);
                        end
                    end
                end
                valueEnter = par.beta.*squeeze(v(:,1,1,:,:,y+1));
                policyEnter(:,:,:,y) = exp(valueEnter./par.wrk.rho)./sum(sum(exp(valueEnter./par.wrk.rho),3),2);
            end
        end
        
        %% Solve for path of wages that clear excess labor demand functions
        function [wages1, frmDensity1, wrkDensity1, errWages, errDensityFrm, errsupplyLabor] = geq(env,par,init,sol,cRobot,wages0,frmDensity0,wrkDensity0,tStart)
            
            %% Equilibrium variables
            % Cost grid and prices with productivity adjustments + Manufacturing demand for profit grid
            frmCost = permute(repmat(par.costAdj(1,:),par.frm.zPoints,1,2),[1 3 2]).*func.cost(env,par,squeeze(wages0(:,1,:)));
            frmPrice = func.priceIndex(par,frmCost,frmDensity0);
            
            % Calculate manufacturing demand for into profit grid
            supplyLabor0 = permute(repmat(init.wrk.mass,env.nOcc,1,env.nSectors),[1 3 2]).*squeeze(sum(sum(sum(wrkDensity0.*par.wrk.hcap,3),2),1));
            demand0 = func.demand(env,par,frmDensity0,supplyLabor0,wages0);
            profit0 = func.profitGrid(par,frmCost,frmPrice,demand0(1,:));
            
            % Consumer prices
            cpi0 = func.cpi(env,par,wages0,frmDensity0);
            
            %% Firms
            % Solve firm's dynamic program
            frmV = solve.frm(env,par,sol,profit0./permute(repmat(cpi0,2,1,par.frm.zPoints), [3 1 2]),cRobot./cpi0,tStart);
            
            % Simulate firm densities
            [frmDensity1, ~] = simulate.frm(env,par,frmDensity0,cRobot./cpi0,frmV,tStart);
            
            %% Workers
            % Worker dynamic programming
            [~, wrkPolicy, wrkPolicyEnter] = solve.wrk(env,par,wages0./permute(repmat(cpi0,env.nOcc+1,1,env.nSectors),[1 3 2]),tStart); % wrk(env,par,wages,tStart)
            
            % Simulate worker densities
            wrkDensity1= simulate.wrk(env,init,wrkDensity0,wrkPolicy,wrkPolicyEnter,tStart);
            
            %% Equilibrium
            supplyLabor1 = permute(repmat(init.wrk.mass,env.nOcc,1,env.nSectors),[1 3 2]).*squeeze(sum(sum(sum(wrkDensity1.*par.wrk.hcap,3),2),1));
            demand1 = func.demand(env,par,frmDensity1,supplyLabor0,wages0);
            factorShareManuf = func.AggFactorShare(env,par,frmDensity1,squeeze(wages0(:,1,:)));
            
            wbAgg = nan(env.nOcc,env.nSectors,env.nYears);
            wbAgg(:,2,:) = squeeze(par.alpha(1:end-1,2,:)).*demand1(2,:);
            wbAgg(:,1,:) = (factorShareManuf(1:end-1,:).*demand1(1,:))./par.frm.markup;
            
            wages1 = wages0;
            wages1(1:end-1,:,tStart+1:end) = wages0(1,2,tStart+1:end).*((wbAgg(:,:,tStart+1:end)./wbAgg(1,2,tStart+1:end))./(supplyLabor1(:,:,tStart+1:end)./supplyLabor1(1,2,tStart+1:end))); % Normalize by occ1-sec2
            errWages = sum(sum(sum(abs(wages1(1:end-1,:,:)-wages0(1:end-1,:,:)))))
            errDensityFrm = sum(sum(sum(abs(frmDensity1-frmDensity0))))
            errsupplyLabor = sum(sum(sum(abs(supplyLabor1-supplyLabor0))))
        end
        
        %% Solve for path of wages that clear excess labor demand functions
        function [wages1, frmDensity1, errWages, errDensityFrm] = geqExL(supplyLabor,env,par,sol,cRobot,wages0,frmDensity0,tStart)
            
            %% Equilibrium variables
            % Cost grid and prices with productivity adjustments + Manufacturing demand for profit grid
            frmCost = permute(repmat(par.costAdj(1,:),par.frm.zPoints,1,2),[1 3 2]).*func.cost(env,par,squeeze(wages0(:,1,:)));
            frmPrice = func.priceIndex(par,frmCost,frmDensity0);
            
            % Calculate manufacturing demand for into profit grid
            demand0 = func.demand(env,par,frmDensity0,supplyLabor,wages0);
            profit0 = func.profitGrid(par,frmCost,frmPrice,demand0(1,:));
            
            % Consumer prices
            cpi0 = func.cpi(env,par,wages0,frmDensity0);
            
            %% Firms
            % Solve firm's dynamic program
            frmV = solve.frm(env,par,sol,profit0./permute(repmat(cpi0,2,1,par.frm.zPoints), [3 1 2]),cRobot./cpi0,tStart);
            
            % Simulate firm densities
            [frmDensity1, ~] = simulate.frm(env,par,frmDensity0,cRobot./cpi0,frmV,tStart);
            
          
            %% Equilibrium
            demand1 = func.demand(env,par,frmDensity1,supplyLabor,wages0);
            factorShareManuf = func.AggFactorShare(env,par,frmDensity1,squeeze(wages0(:,1,:)));
            
            wbAgg = nan(env.nOcc,env.nSectors,env.nYears);
            wbAgg(:,2,:) = squeeze(par.alpha(1:end-1,2,:)).*demand1(2,:);
            wbAgg(:,1,:) = (factorShareManuf(1:end-1,:).*demand1(1,:))./par.frm.markup;
            
            wages1 = wages0;
            wages1(1:end-1,:,tStart+1:end) = wages0(1,2,tStart+1:end).*((wbAgg(:,:,tStart+1:end)./wbAgg(1,2,tStart+1:end))./(supplyLabor(:,:,tStart+1:end)./supplyLabor(1,2,tStart+1:end))); % Normalize by occ1-sec2
            errWages = sum(sum(sum(abs(wages1(1:end-1,:,:)-wages0(1:end-1,:,:)))))
            errDensityFrm = sum(sum(sum(abs(frmDensity1-frmDensity0))))
        end     
        
        
    end
end
