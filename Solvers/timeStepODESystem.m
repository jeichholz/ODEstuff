function soln=timeStepODESystem(ODESystem,vars,InitialConditions,tspan,PostProcessFunctions)

    if class(ODESystem)=="symfun"
        ODESystem=ODESystem(t);
    end

    if class(vars)=="symfun"
        vars=vars("t");
    end

   ics=InitialConditions;
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


    [firstOrderSystem,S]=odeToVectorField(ODESystem);

    symToIdx=dictionary();
    %Convert between naming conventions :(
    for i=1:numel(vars)
        v=vars(i);
        vsym=funsToSyms(v);
        symToIdx(string(vsym))=find(S==vsym);
        differentiationOrder=1;
        while (1)
            vv=diff(v,differentiationOrder);
            vvsym=funsToSyms(vv);
            if differentiationOrder==1
                otherName=sym(['D',char(vsym)]);
            elseif differentiationOrder>1
                otherName=sym(['D',num2str(differentiationOrder),char(vsym)]);
            end
            idx=find(S==otherName);
            if isempty(idx)
                break;
            else
                symToIdx(string(vvsym))=idx;
                differentiationOrder=differentiationOrder+1;
            end
        end
    end
    
    Y0=zeros(numel(ics),1);
    %Set initial conditions. 
    for field=reshape(fieldnames(ICstruct),1,[])
        Y0(symToIdx(field{1}))=getfield(ICstruct,field{1});
    end

    %Now solve. 
    [t,ys]=ode23s(matlabFunction(firstOrderSystem,'vars',{'t','Y'}),tspan,Y0);

    %Put it in to a structure. 
    soln=struct();
    soln.t=t;
    for field=reshape(fieldnames(ICstruct),1,[])
        field=field{1};
        soln=setfield(soln,field,ys(:,symToIdx(field)));
    end

    ICstructnames=fieldnames(ICstruct);
    for i=1:numel(PostProcessFunctions)
        LH=funsToSyms(lhs(PostProcessFunctions(i)));
        RH=matlabFunction(funsToSyms(rhs(PostProcessFunctions(i))),'vars',ICstructnames);
        val=0;
        evalstr=['RH(soln.',ICstructnames{1}];
        for j=2:numel(ICstructnames)
            evalstr=[evalstr,',soln.',ICstructnames{j}];
        end
        evalstr=[evalstr,');'];
        val=eval(evalstr);
        soln=setfield(soln,string(LH),val);
    end
end