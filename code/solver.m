%% Create arrays for endonous variables of GE shooting algorithm
geq = struct();

% Initial values of GE variables (initialization method: Start from baseline estimation)
geq.wages0 = permute(repmat(init.wages,1,1,1,2),[4 1 2 3]);
geq.wages1 = geq.wages0;

geq.frm.density0 = permute(repmat(init.frm.density,1,1,1,2),[4 1 2 3]);
geq.frm.density1 = nan(2,par.frm.zPoints,2,env.nYears);

geq.wrk.density0 = permute(repmat(init.wrk.density,1,1,1,1,1,1,2),[7 1 2 3 4 5 6]);
geq.wrk.density1 = zeros(2,env.wrk.nSkills,env.wrk.nAge,env.wrk.nTen,env.nOcc,env.nSectors,env.nYears);

%% Start fixed-point shooting algorithm
exper.title(sp.exper)

errDensityFrm = 1;
errWages = 1;
errsupplyLabor = 1;
iterShooting = 0;
timerGE = tic;

% Shooting 1: Baseline (costNo=1)
costNo = 1

while (errDensityFrm>sol.tol.ge.densityFrm || errsupplyLabor>sol.tol.ge.supplyLabor || errWages>sol.tol.ge.wages) && iterShooting<sol.iter.ge
    iterShooting = iterShooting+1
    timerIter = tic;
    
    [geq.wages1(costNo,:,:,:), geq.frm.density1(costNo,:,:,:), geq.wrk.density1(costNo,:,:,:,:,:,:), errWages, errDensityFrm, errsupplyLabor] = solve.geq(env,par,init,sol,squeeze(exper.cRobot(sp.exper,costNo,:)).',squeeze(geq.wages0(costNo,:,:,:)),squeeze(geq.frm.density0(costNo,:,:,:)),squeeze(geq.wrk.density0(costNo,:,:,:,:,:,:)),1);
    
    geq.wages0(costNo,1:end-1,:,:) = sol.lambda.wages.*geq.wages0(costNo,1:end-1,:,:) + (1-sol.lambda.wages).*geq.wages1(costNo,1:end-1,:,:);
    geq.frm.density0(costNo,:,:,:) = sol.lambda.frm.*geq.frm.density0(costNo,:,:,:) + (1-sol.lambda.frm).*geq.frm.density1(costNo,:,:,:);
    geq.wrk.density0(costNo,:,:,:,:,:,:) = sol.lambda.wrk.*geq.wrk.density0(costNo,:,:,:,:,:,:) + (1-sol.lambda.wrk).*geq.wrk.density1(costNo,:,:,:,:,:,:);
        
    toc(timerIter)
end

% Shooting 2: Surprise (costNo=2)
if exper.surprise(sp.exper)==1
    costNo = 2
    
    % Baseline equilibrium as initial condition and guess
    geq.frm.density0(costNo,:,:,:) = geq.frm.density0(costNo-1,:,:,:);
    geq.wrk.density0(costNo,:,:,:,:,:,:) = geq.wrk.density0(costNo-1,:,:,:,:,:,:);
    geq.wages0(costNo,:,:,:) = geq.wages0(costNo-1,:,:,:);
    if sp.exper==1
        geq.frm.density0(costNo,:,:,:) = repmat(init.frm.density(:,:,1),1,1,length(init.frm.density)); % Set "zero adoption" as initial guess for "No Robots" counterfactual
    end
    
    errDensityFrm = 1;
    errWages = 1;
    errsupplyLabor = 1;
    iterShooting = 0;
    timerGE = tic;
    
    while (errDensityFrm>sol.tol.ge.densityFrm || errsupplyLabor>sol.tol.ge.supplyLabor || errWages>sol.tol.ge.wages) && iterShooting<sol.iter.ge
        iterShooting = iterShooting+1
        timerIter = tic;
        
        [geq.wages1(costNo,:,:,:), geq.frm.density1(costNo,:,:,:), geq.wrk.density1(costNo,:,:,:,:,:,:), errWages, errDensityFrm, errsupplyLabor] = solve.geq(env,par,init,sol,squeeze(exper.cRobot(sp.exper,costNo,:)).',squeeze(geq.wages0(costNo,:,:,:)),squeeze(geq.frm.density0(costNo,:,:,:)),squeeze(geq.wrk.density0(costNo,:,:,:,:,:,:)),exper.tSurprise(sp.exper));
        
        geq.wages0(costNo,1:end-1,:) = sol.lambda.wages.*geq.wages0(costNo,1:end-1,:) + (1-sol.lambda.wages).*geq.wages1(costNo,1:end-1,:);
        geq.frm.density0(costNo,:,:,:) = sol.lambda.frm.*geq.frm.density0(costNo,:,:,:) + (1-sol.lambda.frm).*geq.frm.density1(costNo,:,:,:);
        geq.wrk.density0(costNo,:,:,:,:,:,:) = sol.lambda.wrk.*geq.wrk.density0(costNo,:,:,:,:,:,:) + (1-sol.lambda.wrk).*geq.wrk.density1(costNo,:,:,:,:,:,:);
        
        toc(timerIter)
    end
    toc(timerGE);
end

%% Save output
geq.cRobot = squeeze(exper.cRobot(sp.exper,:,:));
geq.wages = geq.wages0;
geq.frm.density = geq.frm.density0;
geq.wrk.density = geq.wrk.density0;
geq = rmfield(geq,{'wages0','wages1'});
geq.frm = rmfield(geq.frm,{'density0','density1'});
geq.wrk = rmfield(geq.wrk,{'density0','density1'});
%%
geq = output.geVar(env,par,init,sol,geq);
save(char(strcat('../output/geq',filename,exper.title(sp.exper),'.mat')), '-struct', 'geq');