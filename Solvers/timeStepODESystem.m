function soln=timeStepODESystem(ODESystem,tspan,InitialConditions,PostProcessFunctions,doPlot,doDiagnostics)
%function soln=timeStepODESystem(ODESystem,vars,InitialConditions,tspan,PostProcessFunctions)
%
%Take the system of ODE defined in ODESystem and plug it in to ODE23s.
%Optionally, evaluate the quantities of interest defined in
%PostProcessFunctions. 
%
%INPUT
%
%ODESystem -- a system of pure ODEs, no DAEs.  Should be given as a vector
%of symbolic functions/expressions.  Not a cell. 
%
%InitialConditions -- A cell array of initial conditions.  Initial conditions can be
%specified in a number of different ways.  Examples:
%   x==1  <-- set symbolic function x to a value
%   x(t) == 1 <-- set expression x(t) to a value
%   diff(x) == 1 <-- set derivative of symbolic function x to 1. 
%   "Dxt==1" <-- set derivative of x to 1.  The naming convention is
%   "D<function name><derivative order>.  So "Dthetatt==2" means the second
%   derivative of theta with respect to t. 
%
%tspan -- a vector containing the time span to integrate over. 
%
%PostProcessingFunctions -- optional.  If none, omit or set to [].  If
%included, a list of explicit equations detailing how to calculate
%quantities of interest from the state functions.  For instance, in a
%single pendulum, if the state variables are theta and diff(theta), one
%post-processing function might be 
%x==sin(theta). 
%
%doPlot -- optional.  If omitted or set to [] then defaults to 1.  If the
%value is 1 then a plot of all calculated quantities is calculated. 
%
%%doDiagnositics -- Optional, with default value 1.  If true, then print
%diagnostic information as you proceed. 
%
%OUTPUT
%
%soln -- a struct with one field for each quantity of interest. If
%postprocessing equations are included, those fields are included in soln. 

    if ~exist("doPlot","var") || isempty(doPlot)
        doPlot=1;
    end

    if ~exist("doDiagnostics","var") || isempty(doDiagnostics)
        doDiagnostics=1;
    end

    ODESystem=evalSymfun2Expr(ODESystem);
    InitialConditions=evalSymfun2Expr(InitialConditions);

    ICstruct=parseConditions(InitialConditions);

    [firstOrderSystem,S]=odeToVectorField(ODESystem);

    [~,~,v2f,~,baseFunctions]=symFunsToSymVars(ODESystem);

    vars=unique(v2f(baseFunctions.values));

    dfprintf("The unknown functions (just the base functions, not derivatives) are: %s",vars(1));
    for i=2:numel(vars)
        dfprintf(",%s",vars(i));
    end
    dfprintf("\n")

    %Create a map of state variables to indices. 
    varToIdx=createVarToIndexMap(vars,S);
    
    %Create initial condition vector. 
    Y0=structToVec(ICstruct,varToIdx);

    %Now solve. 
    [t,ys]=ode23s(matlabFunction(firstOrderSystem,'vars',{'t','Y'}),tspan,Y0);

    %Put it in to a structure. 
    soln=struct();
    soln.t=t;
    for field=reshape(fieldnames(ICstruct),1,[])
        field=field{1};
        soln=setfield(soln,field,ys(:,varToIdx(field)));
    end

    %Massivly ugly way to evaluate all the post-processing functions 
    %using the data in the solved ODEs.  
    ICstructnames=fieldnames(ICstruct);


    %Build a string to evaluate the function RH with the fields in soln
    %plugged in. 
    evalRHString=['RH(soln.t,soln.',ICstructnames{1}];
    for j=2:numel(ICstructnames)
        evalRHString=[evalRHString,',soln.',ICstructnames{j}];
    end
    evalRHString=[evalRHString,');'];

    %Create a new RH for each different PP equation, and evaluate. 
    if exist("PostProcessFunctions","var") && ~isempty(PostProcessFunctions)
        PP=symFunsToSymVars(PostProcessFunctions);
        for i=1:numel(PP)
            LH=lhs(PP(i));
            RH=matlabFunction(rhs(PP(i)),'vars',{'t',ICstructnames{:}});
            val=eval(evalRHString);
            soln=setfield(soln,string(LH),val);
        end
    end

    %Plotting!
    if doPlot
        %Plot yo-self!
        figure();
        hold on
        for i=1:numel(ICstructnames)
            plot(soln.t,getfield(soln,ICstructnames{i}),'LineWidth',2)
        end
        legend(ICstructnames,'Location','Best')

        grid on
        hold off
    end


    function dfprintf(varargin)
        if doDiagnostics
            fprintf(varargin{:})
        end
    end


    function symToIdx=createVarToIndexMap(vars,S)
        [~,~,varsToFuns,~,baseVars]=symFunsToSymVars(vars);
        baseVars=unique(baseVars.values);
        symToIdx=dictionary();
        %Convert between naming conventions :(
        for ii=1:numel(baseVars)
            v=baseVars(ii);
            symToIdx(string(v))=find(S==v);
            differentiationOrder=1;
            while (1)
                vv=sym(['D',char(v),repmat('t',1,differentiationOrder)]);
                if differentiationOrder==1
                    otherName=sym(['D',char(v)]);
                elseif differentiationOrder>1
                    otherName=sym(['D',num2str(differentiationOrder),char(v)]);
                end
                idx=find(S==otherName);
                if isempty(idx)
                    break;
                else
                    symToIdx(string(vv))=idx;
                    differentiationOrder=differentiationOrder+1;
                end
            end
        end
    end



end