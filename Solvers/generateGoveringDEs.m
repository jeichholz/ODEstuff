function [systemOfODEs,postProcessingEqns]=generateGoveringDEs(DEs,positionConstraints,velocityConstraints,...
    stateVariables,parameters,doDiagnostics)
%function [systemOfODEs,postProcessingEqns]=analyzeSystem(DEs,positionConstraints,velocityConstraints, ...
%    unknownFunctions,stateVariables)
%
%Take a system of DAE describing rigid-body motion, and try to find the
%governing system of pure ODE, along with a list of "post-processing"
%functions which describe all other quantities of interest if the state
%functions are known.  
%
%INPUTS
%
%DEs -- A system of DEs, and *other equations that don't need to be
%differentiated*.  This would generally include mx''=F and Ig theta''=M
%type equations, as well as other constituatitive relationns like F1=-F3.
%This should be a vector and not, e.g., a cell. 
%
%PositionConstraints -- optional.  If there are none then either omit or
%set to [].  A list of constraints that should be differentiated twice.
%These would tend to include only position-like functions. 
%
%VelocityConstraints -- optional.  If there are none then either omit or
%set to [].  A list of constraints that need to be differentiated once.
%These would tend to include only velocity-like functions. 
%
%unkownFunctions -- *not optional*.  The unknown time-varying function to solve for. Do not include paramters.
% Only include the "base" functions.  So include theta, but not diff(theta). 
%
%stateVariables -- The state variables that you which to express the system
%in terms of.  This tool can warn you if you give the wrong number of state
%variables. For instance a single pendulum needs two state variables, you will be warned if you supply 1 or 3. 
%However, supplying the right number of state variables is not enough. For
%instance, supplying x(t) and theta(t) as state variables for a single
%pendulum will not suceed, but x(t) and diff(theta) will. Note that if you
%supply the "wrong" state variables, this tool may try to work for a *very*
%long time.  In that case, use ctrl+c and try again. 
%
%parameters -- a vector of symbols which should be treated as paramters and
%not unknown functions.  For the most part, non time-varying variables are correctly treated as paramters by default.  
%However, generic non-homogeneous terms, like the f(t) in x''+x=f(t), need
%to be specified as parameters. 
%
%doDiagnositics -- Optional, with default value 1.  If true, then print
%diagnostic information as you proceed. 
%
%OUTPUTS
%
%systemOfODEs -- a system of pure ODEs in the supplied state variables that
%characterize the original system. 
%
%postProcessingEqns -- a system of equations that tell you how to calculate
%all other quantities of interest once your state functions are
%approximated. 
%
%CAVEATS
%
% *The time variable must be called t. 


    if ~exist("doDiagnostics","var") || isempty(doDiagnostics)
        doDiagnostics=1;
    end

    if ~exist("parameters","var") || isempty(parameters)
        parameters=[];
    end

    %Make sure the inputs are expressions, not functions. Cast them if
    %needed. 
    DEs=evalSymfun2Expr(DEs);
    positionConstraints=evalSymfun2Expr(positionConstraints);
    velocityConstraints=evalSymfun2Expr(velocityConstraints);
    parameters=evalSymfun2Expr(parameters);
    
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
    

    [system,funsToVars,varsToFuns,order,baseFunctions]=symFunsToSymVars(system);

    if ~isempty(parameters)
        parameters=symFunsToSymVars(parameters);
    end
    
    unknownFunctions=setdiff(unique(baseFunctions.values),parameters);

    dfprintf('The unknown functions are:\n')
    for i=1:numel(unknownFunctions)
       dfprintf("%s",string(unknownFunctions(i)));
       if  order(symFunsToSymVars(unknownFunctions(i)))==0
            dfprintf('(force)');
       end

       dfprintf("\n")
    end
    dfprintf("\n");

    numeqns=numel(system);
    numsymbols=numel(symvar(system));
    numparameters=numsymbols-numeqns;

    numstatevarformula=sum(order.values)-numel(velocityConstraints)-2*numel(positionConstraints);
    dfprintf("Expected Number of State Variables:%d\n",numstatevarformula);
    

    %If the user told us what should be the "state variables", use them.
    %Otherwise, try to find them. 
    if exist('stateVariables','var') && ~isempty(stateVariables)
        stateVariables=evalSymfun2Expr(stateVariables);

        if numel(stateVariables) ~= numstatevarformula
            warning('\nWarning.  This system will need to be expressed in terms of %d state variables.  You have provided %d.\n',numstatevarformula,numel(stateVariables));
        end

        soln=solveInTermsOf(system,symFunsToSymVars(stateVariables));
        
        %Which equations to grab?  All of the ones that are the highest
        %derivative of a state function. 
        solveInTermsOfBase=unique(baseFunctions(symFunsToSymVars(stateVariables)));
        
        allLHS=sym(fieldnames(soln));
        systemOfODELHS=sym([]);
        systemOfODEs=sym([]);
        postProcessingEqns=sym([]);

        for i=1:numel(solveInTermsOfBase)
            f=solveInTermsOfBase(i);
            LHS=funsToVars(diff(varsToFuns(f),order(f)));
            systemOfODELHS(end+1)=LHS;
        end
        
        postProcessingEqnsLHS=setdiff(allLHS,systemOfODELHS);
        
        %The LHS (second derivatives) are now in systemOfODESLHS and
        %postProcessingLHS.  Now just build the equations. 
        for i=1:numel(systemOfODELHS)
            if numel(systemOfODELHS)<=2
                systemOfODEs(end+1,1)=systemOfODELHS(i)==simplify(selectBranchOrCrash(soln,systemOfODELHS(i)));
            else
                systemOfODEs(end+1,1)=systemOfODELHS(i)==selectBranchOrCrash(soln,systemOfODELHS(i));
            end
        end

        for i=1:numel(postProcessingEqnsLHS)
            postProcessingEqns(end+1,1)=postProcessingEqnsLHS(i)==selectBranchOrCrash(soln,postProcessingEqnsLHS(i));
        end


        %Now reconstitute back to symfuns, not symvars.
        systemOfODEs=subs(systemOfODEs,varsToFuns.keys,varsToFuns.values);
        postProcessingEqns=subs(postProcessingEqns,varsToFuns.keys,varsToFuns.values);


    else
        %In this case, we'll loop over all "reasonable" choices of state
        %variables.  
        error("Determining state variables without a hint is not currently supported.\n")
        %WARNING: THIS MIGHT BE SLOW IF THERE ARE LOTS OF UNKOWN FUNCTIONS.
        for i=1:numel(unknownDifferentialFunctions)
            possibleCombos=nchoosek(unknownDifferentialFunctions,i);
            for j=1:size(possibleCombos,1)
                totalOrder=0;
                for k=1:i
                    totalOrder=totalOrder+order(possibleCombos(j,k));
                end
                if totalOrder==numparameters
                    dfprintf("It seems resonable to express everything in terms of: ")
                    stateVars=sym([]);
                    for k=1:i
                        dfprintf(string(possibleCombos(j,k)));
                        stateVars(end+1)=possibleCombos(j,k);
                        for l=1:order(possibleCombos(j,k))-1
                            stateVars(end+1)=diff(possibleCombos(j,k),l);
                        end
                    end
                    dfprintf("\n")
                    soln=solveInTermsOf(system,[reshape(stateVars,1,[]),reshape(paremeters,1,[])]);

                end
            end
        end
    end


    function dfprintf(varargin)
        if doDiagnostics
            fprintf(varargin{:});
        end
    end

    function RHS=selectBranchOrCrash(solnStruct,field)
            RHS=getfield(solnStruct,string(field));
            if isempty(RHS)
                warning("Error! Couldn't solve for: %s.  Probably going to crash now.\n",string(field));
            end
            if numel(RHS)>1
                dfprintf("Warning! There are multiple branches for %s. This likely means that this dynamical system is piecewise-defined in terms of the state variables.  Arbitrarily selecting the first branch.\n",string(field));
                RHS=RHS(1);
            end
    end

    function soln=solveInTermsOf(system,stateFunctions)
        solveFor=setdiff(setdiff(varsToFuns.keys,stateFunctions),parameters);
        soln=solve(system,solveFor);
    end


end