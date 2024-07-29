function [systemOfODEs,postProcessingEqns]=analyzeSystem(DEs,positionConstraints,velocityConstraints, ...
    unknownFunctions,stateVariables)
    
    %Make sure the inputs are expressions, not functions. Cast them if
    %needed. 
    DEs=evalSymfun2Expr(DEs);
    positionConstraints=evalSymfun2Expr(positionConstraints);
    velocityConstraints=evalSymfun2Expr(velocityConstraints);
    unknownFunctions=evalSymfun2Expr(unknownFunctions);
    
    %Position constraints are constraints that need to be differentiated
    %twice. Add all position constraints to the system. 
    system=DEs;
    for i=1:numel(positionConstraints)
        system(end+1)=positionConstraints(i);
        system(end+1)=diff(positionConstraints(i));
        system(end+1)=diff(positionConstraints(i),2);
    end

    %Velocity constraints are constraints that need to be differentiated
    %once.  Add all position constraints to the system. 
    for i=1:numel(velocityConstraints)
        system(end+1)=velocityConstraints(i);
        system(end+1)=diff(velocityConstraints(i));
    end
    
    fprintf('There are %d equations, and %d unknown time-varying functions.\n',numel(system),numel(unknownFunctions));

    order=dictionary();

    for i =1:numel(unknownFunctions)
        order(unknownFunctions(i))=0;
        for j=1:5
            if any(has(system,diff(unknownFunctions(i),j)))
                order(unknownFunctions(i))=j;
            end
        end
    end

    [system,funcs2str,str2funcs]=funsToSyms(system);

    unknownDifferentialFunctions=sym([]);
    fprintf('The algebraic (force) functions are:\n')
    for i=1:numel(unknownFunctions)
       if order(unknownFunctions(i))>=1
            unknownDifferentialFunctions(end+1)=unknownFunctions(i);
       else
            fprintf('%s\n',string(unknownFunctions(i)));
       end
    end
    fprintf("\n");

    numeqns=numel(system);
    numsymbols=numel(symvar(system));
    numparameters=numsymbols-numeqns;
    fprintf("Adding in all constraints, and counting pure symbols (so that x, Dxt, Dxtt count as 3 symbols) I get:\n")
    fprintf("equations: %d\nsymbols:%d\nfree parameters: %d.\n\n", ...
        numel(system),numel(symvar(system)),numparameters)
    fprintf("Expected number of state variables: %d\n\n",numparameters)

    %If the user told us what should be the "state variables", use them.
    %Otherwise, try to find them. 
    if exist('stateVariables','var') && ~isempty(stateVariables)
        if class(stateVariables)=="symfun"
            stateVariables=stateVariables("t");
        end
        if numel(stateVariables) ~= numparameters
            fprintf('Warning.  This system will need to be expressed in terms of %d paramters.  You have provided %d.\n',numparameters,numel(stateVariables));
        end
        soln=solveInTermsOf(system,stateVariables);
        %Ok, now split it off into the pure ODE part and then the
        %post-processing part. 
        targetList=sym([]);
        for i=1:numel(stateVariables)
            %differentiate the state variables until you get something that
            %you solved for.  
            target=stateVariables(i);
            while any(ismember(stateVariables,target))
                target=diff(target);
            end
            targetList(end+1)=target;
        end
        targetList=unique(targetList);
        postList=fieldnames(soln);
        for i=1:numel(targetList)
            postList=setdiff(postList,funcs2str(targetList(i)));
        end
        systemOfODEs=sym([]);
        postProcessingEqns=sym([]);
        for i=1:numel(targetList)
            RHS=getfield(soln,funcs2str(targetList(i)));
            if isempty(RHS)
                fprintf("Error! Couldn't solve for: %s.  Probably going to crash now.\n",string(targetList(i)));
            end
            if numel(RHS)>1
                fprintf("Warning! There are multiple branches for %s. This likely means that this dynamical system is piecewise-defined in terms of the state variables.  Arbitrarily selecting the first branch.\n",string(targetList(i)));
                RHS=RHS(1);
            end
            systemOfODEs(end+1,1)=funcs2str(targetList(i))==simplify(RHS);
        end
        for i=1:numel(postList)
            RHS=getfield(soln,postList(i));
            if isempty(RHS)
                fprintf("Error! Couldn't solve for: %s.  Probably going to crash now.\n",string(postList(i)));
            end
            if numel(RHS)>1
                fprintf("Warning! There are multiple branches for %s. This likely means that this dynamical system is piecewise-defined in terms of the state variables.  Arbitrarily selecting the first branch.\n",string(postList(i)));
                RHS=RHS(1);
            end
            postProcessingEqns(end+1,1)=sym(postList(i))==simplify(RHS);
        end
        for f=reshape(unknownFunctions,1,[])
            systemOfODEs=subs(systemOfODEs,funcs2str(f),f);
            postProcessingEqns=subs(postProcessingEqns,funcs2str(f),f);
            for i=1:order(f)
                systemOfODEs=subs(systemOfODEs,funcs2str(diff(f,i)),diff(f,i));
                postProcessingEqns=subs(postProcessingEqns,funcs2str(diff(f,i)),diff(f,i));            
            end
        end

    else
        %In this case, we'll loop over all "reasonable" choices of state
        %variables.  

        %WARNING: THIS MIGHT BE SLOW IF THERE ARE LOTS OF UNKOWN FUNCTIONS.
        for i=1:numel(unknownDifferentialFunctions)
            possibleCombos=nchoosek(unknownDifferentialFunctions,i);
            for j=1:size(possibleCombos,1)
                totalOrder=0;
                for k=1:i
                    totalOrder=totalOrder+order(possibleCombos(j,k));
                end
                if totalOrder==numparameters
                    fprintf("It seems resonable to express everything in terms of: ")
                    stateVars=sym([]);
                    for k=1:i
                        fprintf(string(possibleCombos(j,k)));
                        stateVars(end+1)=possibleCombos(j,k);
                        for l=1:order(possibleCombos(j,k))-1
                            stateVars(end+1)=diff(possibleCombos(j,k),l);
                        end
                    end
                    fprintf("\n")
                    soln=solveInTermsOf(system,stateVars);

                end
            end
        end
    end


    function soln=solveInTermsOf(system,stateFunctions)
        tmp=strings(numel(stateFunctions),1);
        for ii=1:numel(stateFunctions)
            tmp(ii)=funcs2str(stateFunctions(ii));
        end
        solveFor=setdiff(str2funcs.keys,tmp);
        soln=solve(system,cellstr(solveFor));
    end


end