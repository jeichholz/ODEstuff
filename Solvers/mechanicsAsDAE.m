function [soln]=mechanicsAsDAE(DEs,PositionConstraints,VelocityConstraints,tspan,ics,hints,doPlot,doDiagnostics,solverChoice)
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
%solverChoice -- Optional.  Set to one of "ode15i", "ode23t", or "ode15s".
%Default value is ode15s. 


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

    %If the DEs are a symfun, just turn them into an expression. 
    if class(DEs)=="symfun"
        DEs=DEs("t");
    end

    if class(PositionConstraints)=="symfun"
        PositionConstraints=PositionConstraints("t");
    end
    
    if ~exist("solverChoice","var")||isempty(solverChoice)
        solverChoice="ode15s";
    else
        validSolverChoices={'ode15s','ode23t','ode15i'};
        if ~ismember(solverChoice,validSolverChoices)
            error('Solver choice, %s, must be one of: %s.\n',solverChoice,sprintf("%s ",validSolverChoices{:}));
        end
    end


    if class(VelocityConstraints)=="symfun"
        VelocityConstraints=VelocityConstraints("t");
    end

    if class(vars)=="symfun"
        vars=vars("t");
    end

    %Sometimes the initial conditions might be interpreted as a
    %vector-valued symfun.  Stop that. 
    if class(ics)=="symfun"
        ics=ics("t");
    end

   
    %It is reasonable to pass in initial conditions as a vector of symbols,
    %but not what we want.  Convert to cell array. 
    if class(ics)=="sym"
        ics=sym2cell(ics);
    end

    %Store the initial conditions in a structure for later
    %substitution/manipulation. 
    ICstruct=struct();

    %Initial conditions may be in the form of a string, "y==4", "Dyt==4"
    %Character array,                                   'y==4', 'Dyt==4'
    %symbolic function value,                           y==4,
    %diff(y,t)==4
    %or symvol value,                                   y(t)==4, diff(y(t),t)==4 
    for i=1:numel(ics)
        if class(ics{i})=="string"
            L=ics{i}.extractBefore("==");
            R=sym(ics{i}.extractAfter("=="));
        elseif class(ics{i})=="char"
            equalSigns=find(ics{i}=="=");
            L=ics{i}(1:equalSigns(1)-1);
            R=sym(ics{i}(equalSigns(2)+1:end));
        else
            [~,asSymbol]=funsToSyms(ics{i});
            if class(lhs(ics{i}))=="symfun"
                tmp=lhs(ics{i});
                L=asSymbol(tmp("t"));
                R=rhs(ics{i});
                R=R("t");
            else
                L=asSymbol(lhs(ics{i}));
                R=rhs(ics{i});
            end
        end
        ICstruct=setfield(ICstruct,L,R);
    end


    %The list of DE's to solve.  We'll differentiate the constraints
    %the correct number of times, and append them on the to the end of the
    %list. 
    DEs=DEs;

    %Turn the position constraints into two constraintes each, and the
    %velocity constraintes into one each; 
    constraints=sym([]);
    for i=1:numel(PositionConstraints)
        constraints(end+1:end+2,1)=funsToSyms([PositionConstraints(i); diff(PositionConstraints(i),"t")]);
        DEs(end+1)=diff(PositionConstraints(i),"t","t");
    end

    for i=1:numel(VelocityConstraints)
        constraints(end+1,1)=funsToSyms(VelocityConstraints(i));
        DEs(end+1)=diff(VelocityConstraints(i),"t");
    end


    %Reduce it to 1st order. 
    [firstOrder,firstOrderVars] = reduceDifferentialOrder(DEs,vars);

    [~,var2Name,name2Var]=funsToSyms(firstOrderVars);
    varName2Idx=dictionary();
    var2Idx=dictionary();
    idx2VarName={};
    for i = 1:numel(firstOrderVars)
        symbol=firstOrderVars(i);
        varName2Idx(var2Name(symbol))=i;
        var2Idx(symbol)=i;
        idx2VarName{i}=var2Name(symbol);
    end



    %Are there independent blocks?
    %findDecoupledBlocks(firstOrder,firstOrderVars)

    %You need to find initial conditions that satisfy the constraints that
    %were spit out of reduceDAEtoODE. 

    specificConstraints=constraints;

    %ASSUMES THE TIME VARIABLE IS NAMED t
    specificConstraints=subs(specificConstraints,"t",t0);

    %Loop through the initial conditions and do the assignment. IDK why we can't just sub in the ICstruct, but that 
    %doesn't seem to work. 
    fields=fieldnames(ICstruct);
    for i=1:numel(fields)
        field=fields{i};
        specificConstraints=subs(specificConstraints,field,getfield(ICstruct,field));
end

    ICsoln=solve(specificConstraints);
    
    %Loop over all fields of the solution to the constraints.
    fields=fieldnames(ICsoln);
    for fidx=1:numel(fields)
        field=fields{fidx};
        if isempty(getfield(ICsoln,field))
            fprintf("Error: Could not find a suitable initial value for %s.\n",field);
        end
        if numel(getfield(ICsoln,field))>1
            possVals=getfield(ICsoln,field);
            fprintf("Found multiple suitable values for %s:\n ",field)
            for j=1:numel(getfield(ICsoln,field))
                fprintf("\t%f\n",possVals(j));
            end
            fprintf("Arbitrarily selecting the first one.\n")
            ICsoln=setfield(ICsoln,field,double(possVals(1)));
end
end

    %Now we have satisfied the initial constraints of the system, however,
    %that might not have pinned down all of the necessary initial
    %conditions, like unknown forces which will only show up in the DEs
    %themselves. First, build a partial initial condition vector out of the
    %ICs that we were given and that we determined in order to be feasible.
    %Then, use decic to find anything else that we need. 
    Y0=zeros(numel(firstOrderVars),1);
    Y0Fixed=zeros(size(Y0));
    
    fields=fieldnames(ICstruct);
    for i=1:numel(fields)
        field=string(fields{i});
        idx=varName2Idx(field);
        Y0(idx)=getfield(ICstruct,field );
        Y0Fixed(idx)=1;
end
    fields=fieldnames(ICsoln);
    for i=1:numel(fields)
        field=string(fields{i});
        idx=varName2Idx(field);
        Y0(idx)=getfield(ICsoln,field);
        Y0Fixed(idx)=1;
    end

    [M,F] = massMatrixForm(firstOrder,firstOrderVars);
    M = odeFunction(M, firstOrderVars);
    F = odeFunction(F, firstOrderVars);
    implicitForm=@(t,y,yp) M(t,y)*yp-F(t,y);

    [Y0,YP0]=decic(implicitForm,t0,Y0,Y0Fixed,zeros(size(Y0)),zeros(size(Y0Fixed)));


    fprintf("Initial Conditions:\n")
    fprintf("Y0/YP0:\n")
    for i=1:numel(Y0)
        fprintf("%s: %f %s:%f\n",idx2VarName{i},Y0(i),[char(idx2VarName{i}),'t'],YP0(i));
    end

    %Ok, solve!
    opt = odeset('Mass', M,'RelTol', 10.0^(-7), 'AbsTol' , 10.0^(-7),'InitialSlope',YP0);

    if string(solverChoice)=='ode15i'
        dfprintf('Using ode15i solver.\n')
        [tSol,ySol]=ode15i(implicitForm,[t0,tmax],Y0,YP0,opt);
    elseif string(solverChoice)=='ode23t'
        dfprintf('using ode23t solver.\n')
        [tSol,ySol] = ode23t(F, [t0, tmax], Y0, opt);
    elseif string(solverChoice)=='ode15s'
        [tSol,ySol] = ode15s(F, [t0, tmax], Y0, opt);
    end

    soln=struct();
    soln.t=tSol;
    for i=1:numel(idx2VarName)
        soln=setfield(soln,idx2VarName{i},ySol(:,i));
    end

    origVars=numel(vars);

    if ~exist("doPlot","var") || doPlot
        %Plot yo-self!
        plot(tSol,ySol(:,1:origVars),'LineWidth',2)
    
    
        for k = 1:numel(vars)
            S{k} = char(vars(k));
        end
    
        legend(S, 'Location', 'Best')
        grid on
    end