function expr=evalSymfun2Expr(symfun)
    if class(symfun)=="symfun"
        %expr=formula(symfun);
        expr=symfun("t");
    else
        expr=symfun;
    end