function [soln]=mechanicsAsDAE(DEs,PositionConstraints,VelocityConstraints,tspan,ics,hints,doPlot,doDiagnostics)
%function [soln]=mechanicsAsDAE(DEs,PositionConstraints,VelocityConstraints,vars,tspan,ics,hints,doPlot)
%
%Solve a DAE describing some kind of rigid-body motion.  This is purely
%numerical. 
%
%INPUTS: 
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
%tspan -- the time interval to integrate over, in the form [tmin,tmax]. 
%
%ics -- A cell array of initial conditions.  Initial conditions can be
%specified in a number of different ways.  Examples:
%   x==1  <-- set symbolic function x to a value
%   x(t) == 1 <-- set expression x(t) to a value
%   diff(x) == 1 <-- set derivative of symbolic function x to 1. 
%   "Dxt==1" <-- set derivative of x to 1.  The naming convention is
%   "D<function name><derivative order>.  So "Dthetatt==2" means the second
%   derivative of theta with respect to t. 
%
%hints -- Optional, if there are no hints then omit or set to [].  This tool 
%tries to warn you if there are finetely many initial
%states that satisfy your initial conditions.  For example, for a single
%pendulum, if your initial conditions are x==0.5 and diff(x)==0, there are
%two different suitable initial conditions for y.  Use the hints to select
%one of thost initial conditions. The tool will choose the initial
%condition that satisfies the constraints that is closest to your hints.
%Format is identical to ics. 
%
%doPlot -- Optional, if omitted or empty then defaults to 1.  If non-zero then create a plot at the end of the run.  
%
%doDiagnositics -- Options, with default value 1.  If true, then print
%diagnostic information as you proceed. 
%
%
%OUTPUTS
%
%soln -- a structure with a field for every calculated function.
%
%CAVEATS: 
%
% * At the moment the time variable MUST be named t. 
%
% * I intend for it to be ok to pass in either symbolic functions, like x,
% or symbolic expressions like x(t).  If trouble arises, it is safest to
% pass in expressions, like x(t) rather than x, or diff(x(t),t) rather than
% diff(x). 


    t0=tspan(1);
    tmax=tspan(2);

    %Make sure that nothing is a symfun, but rather a vector of symbolic
    %expressions. 
    DEs=evalSymfun2Expr(DEs);
    PositionConstraints=evalSymfun2Expr(PositionConstraints);
    VelocityConstraints=evalSymfun2Expr(VelocityConstraints);

    ICstruct=parseConditions(ics);

    if exist("hints","var") && ~isempty(hints)
        hintsStruct=parseConditions(hints);
    else
        hintsStruct=[];
    end
    
    if ~exist("doPlot","var") || isempty(doPlot)
        doPlot=1;
    end

    if ~exist("doDiagnostics","var") || isempty(doDiagnostics)
        doDiagnostics=1;
    end

    %Turn the position constraints into two constraints each, and the
    %velocity constraints into one each; 
    constraints=sym([]);
    for i=1:numel(PositionConstraints)
        constraints=[constraints;[PositionConstraints(i); diff(PositionConstraints(i),"t")]];
        DEs(end+1)=diff(PositionConstraints(i),"t","t");
    end

    for i=1:numel(VelocityConstraints)
        constraints(end+1,1)=VelocityConstraints(i);
        DEs(end+1)=diff(VelocityConstraints(i),"t");
    end

    [~,~,varsToSyms,~,baseVars]=symFunsToSymVars([DEs;constraints]);

    constraints=symFunsToSymVars(constraints);
    
    vars=unique(varsToSyms(baseVars.values));

    dfprintf("The unknown functions (just the base functions, not derivatives) are: %s",vars(1));
    for i=2:numel(vars)
        dfprintf(",%s",vars(i));
    end
    dfprintf("\n")

    %Reduce the system of DEs to 1st order. 
    [firstOrderSystem,firstOrderFuns] = reduceDifferentialOrder(DEs,vars);

    %The variables in firstOrder are just plain syms, but they have (t) in
    %thier name, which is still an issue for solve.  So rename them to get
    %rid of that. 
    [firstOrderVars,firstOrderFuns2FirstOrderVars]=symFunsToSymVars(firstOrderFuns);
    [~,~,~,diffOrders]=symFunsToSymVars(firstOrderSystem);

    %We'll want to connect those variables to indices in an initial
    %condition soon. 
    firstOrderVarToIdx=dictionary();
    for i = 1:numel(firstOrderVars)
        firstOrderVarToIdx(firstOrderVars(i))=i;
    end

    numConstraints=numel(constraints);
    numSymbols=numel(firstOrderVars);
    numICsGiven=numel(fieldnames(ICstruct));
    numExpectedStates=sum(diffOrders.values)-numConstraints;

    if numICsGiven ~= numExpectedStates
        warning('This system is expected to require %d state variables, but you provided %d initial conditions.\n',numExpectedStates,numICsGiven);
    end

    %You need to find initial conditions that satisfy the constraints that
    %were spit out of reduceDAEtoODE. 

    specificConstraints=constraints;

    %ASSUMES THE TIME VARIABLE IS NAMED t
    specificConstraints=subs(specificConstraints,"t",t0);

    %Substitute the given initial conditions into the constraints. 
    specificConstraints=structsubs(specificConstraints,ICstruct);

    %Solve for what you can symbolically. It is nice to do it symbolically
    %in order to warn the user if there are multiple feasible choices. 
    ICsoln=solve(specificConstraints,'Real',1);
    
    %It is possible that ICsoln has no solutions in some fields, which is
    %bad, or multiple solutions in some fields, in which case we need to
    %look at the hints to decide or decide arbitrarily.  Deal with that. 
    ICsoln=selectIC(ICsoln,hintsStruct);

    %Now we have satisfied the initial constraints of the system, however,
    %that might not have pinned down all of the necessary initial
    %conditions, like unknown forces which will only show up in the DEs
    %themselves. First, build a partial initial condition vector out of the
    %ICs that we were given and that we determined in order to be feasible.
    %Then, use decic to find anything else that we need. 
    
    [Y0,Y0Fixed]=structToVec({ICsoln,ICstruct},firstOrderVarToIdx);

    %Use mass-matrix form for solving the ODE. 
    [Msym,Fsym] = massMatrixForm(firstOrderSystem,firstOrderFuns);

    %For some reason, when the inputs are given as functions like x, rather
    %than x(t), it becomes impossible to use the standard odeFunction to
    %turn M,F into matlab functions.  Convert to symVars rather than
    %anything involving names that even look like functions.  
    M=symFunsToSymVars(Msym);
    M = matlabFunction(M, "Vars",{sym("t"),firstOrderVars});

    F=symFunsToSymVars(Fsym);
    F =matlabFunction(F, "Vars",{sym("t"),firstOrderVars});

    implicitForm=@(t,y,yp) M(t,y)*yp-F(t,y);


    %Use decic to finish out any as-yet determined initial conditions. 
    [Y0,YP0]=decic(implicitForm,t0,Y0,Y0Fixed,zeros(size(Y0)),zeros(size(Y0Fixed)));


    dfprintf("Initial Conditions:\n")
    dfprintf("Y0/YP0:\n")
    for i=1:numel(Y0)
        dfprintf("%s: %f %s:%f\n",string(firstOrderVars(i)),Y0(i),[char(firstOrderVars(i)),'t'],YP0(i));
    end

    %Ok, solve!
    opt = odeset('Mass', M,'RelTol', 10.0^(-7), 'AbsTol' , 10.0^(-7),'InitialSlope',YP0);

    %[tSol,ySol]=ode15i(implicitForm,[t0,tmax],Y0,YP0,opt);
    %[tSol,ySol] = ode23t(F, [t0, tmax], Y0, opt);
    [tSol,ySol] = ode15s(F, [t0, tmax], Y0, opt);


    soln=struct();
    soln.t=tSol;
    for i=1:numel(firstOrderVars)
        soln=setfield(soln,string(firstOrderVars(i)),ySol(:,i));
    end



    if doPlot
        origVars=numel(vars);
        %Plot yo-self!
        plot(tSol,ySol(:,1:origVars),'LineWidth',2)
    
    
        for k = 1:numel(vars)
            S{k} = char(vars(k));
        end
    
        legend(S, 'Location', 'Best')
        grid on
    end



    function dfprintf(varargin)
        if doDiagnostics
            fprintf(varargin{:});
        end
    end


    function val=selectIC(soln,hintsStruct)

        %First, check to see if there are any empty fields and stop. 
        fields=reshape(fieldnames(soln),1,[]);

        for f=fields
            if isempty(getfield(soln,f{1}))
                error("No suitable initial value for %s found.\n",f{1});
            end
        end

        %See if there is a single unique solution. 
        numentries=zeros(numel(fields),1);
        for ii=1:numel(fields)
            numentries(ii)=numel(getfield(soln,fields{ii}));
        end
        
        %If everyone has just one unique solution, great, go with that.
        if all(numentries==1)
            val=soln;
        else
            if ~(all(numentries==numentries(1)))
                error('Error in selecting initial conditions. Possible initial conditions are %s\n',string(soln));
            end
            dfprintf('Multiple feasible initial conditions found:\n')
            for kk=1:numentries(1)
                for jj=1:numel(fields)
                    tmp=getfield(soln,fields{jj});
                    dfprintf('%s: %f    ',fields{jj},tmp(kk));
                end
                dfprintf('\n')
            end

            if ~isempty(hintsStruct)
                %Make an error vector for each possible initial condition. 
                err=zeros(numentries(1),1);
    
                %Loop over each hint entry.  Calculate the distance from each
                %possible IC to the hint. 
                hintFields=fieldnames(hintsStruct);
                for ll=1:numel(hintFields)
                    hf=hintFields{ll};
                    hfv=getfield(hintsStruct,hf);
                    for jj=1:numentries(1)
                        tmp=getfield(soln,hf);
                        err(jj)=err(jj)+(tmp(jj)-hfv)^2;
                    end
                end
    
                %Which possible IC has the lowest error?
                [~,midx]=min(err);

            else
                midx=1;
            end

            dfprintf("\nSelecting IC:\n")
            %Ok, go with midx. 
            for f=fields
                tmp=getfield(soln,f{1});
                soln=setfield(soln,f{1},tmp(midx));
                dfprintf("%s:%f    ",f{1},tmp(midx));
            end
            if isempty(hintsStruct)
                dfprintf('\nBecause it is first.  Use a hint if you want to select a different initial condition.\n\n')
            else
                dfprintf("\nBecause it is closest to satisfying the hints.\n\n")
            end
            val=soln;

                


        end
        
    end


end