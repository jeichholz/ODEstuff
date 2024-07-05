function [soln]=mechanicsAsDAE(DEs,PositionConstraints,VelocityConstraints,vars,t0,tmax,ics,doPlot)

    %If the DEs are a symfun, just turn them into an expression. 
    if class(DEs)=="symfun"
        DEs=DEs("t");
    end

    if class(PositionConstraints)=="symfun"
        PositionConstraints=PositionConstraints("t");
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

    %[tSol,ySol]=ode15i(implicitForm,[t0,tmax],Y0,YP0,opt);
    %[tSol,ySol] = ode23t(F, [t0, tmax], Y0, opt);
    [tSol,ySol] = ode15s(F, [t0, tmax], Y0, opt);


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