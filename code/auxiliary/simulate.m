classdef simulate
    methods(Static)
        
        %% Simulate firm densities (adoption decisions) 
        function [density, adopt0] = frm(env,par,densityInit,cRobot,v,tStart)
            density = nan(par.frm.zPoints,2,env.nYears);
            adopt0 = nan(par.frm.zPoints,env.nYears);
            density(:,:,1:tStart) = densityInit(:,:,1:tStart);
            for t=tStart:1:env.nYears-1
                adopt = (1+exp((cRobot(t)+par.beta.*par.frm.zP*(v(:,1,t+1)-v(:,2,t+1)))/par.frm.nu)).^(-1);
                adopt0(:,t) = adopt;
                density(:,2,t+1) = (1-par.frm.theta).*func.zLom(par,density(:,2,t)) + func.zLom(par,density(:,1,t).*adopt); % NB: we do not take into account that firms are churning through idiosyncratic productivity states via \texttt{par.frm.zP}
                density(:,1,t+1) = func.zLom(par,density(:,1,t).*(1-adopt)) + par.frm.theta.*func.zLom(par,density(:,2,t));
            end
        end
            
        %% Simulate worker densities (occupational choice) 
        function density = wrk(env,init,densityInit,policy,policyEnter,tStart)
            density = zeros(env.wrk.nSkills,env.wrk.nAge,env.wrk.nTen,env.nOcc,env.nSectors,env.nYears);
            density(:,:,:,:,:,1:tStart) = densityInit(:,:,:,:,:,1:tStart);
            
            for y=tStart:env.nYears-1 % Change: tSim --> env.nYears
                % Replace retiring cohort with new cohort
                retiring = sum(sum(sum(sum(sum(sum(density(:,env.wrk.nAge,:,:,:,y)))))));
                density(:,1,1,:,:,y+1) = policyEnter(:,:,:,y).*retiring.*repmat(init.wrk.distEnter(:,y+1),1,env.nOcc,env.nSectors); % Skill distribution of entering cohorts
                for a = 1:env.wrk.nAge-1
                    for t=1:min(a,env.wrk.nTen)
                        for o=1:env.nOcc
                            for s=1:env.nSectors
                                for oNext=1:env.nOcc
                                    for sNext=1:env.nSectors
                                        if oNext==o && sNext==s
                                            density(:,a+1,min(t+1,env.wrk.nTen),oNext,sNext,y+1) = density(:,a+1,min(t+1,env.wrk.nTen),oNext,sNext,y+1) + density(:,a,t,o,s,y).*policy(:,a,t,o,s,oNext,sNext,y);
                                        else
                                            density(:,a+1,1,oNext,sNext,y+1) = density(:,a+1,1,oNext,sNext,y+1) + density(:,a,t,o,s,y).*policy(:,a,t,o,s,oNext,sNext,y);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        
    end
end