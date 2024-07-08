function expr=evalSymfun2Expr(symfun)
    if class(symfun)=="symfun"
        var=symvar(symfun);
        expr=symfun(var);
    else
        expr=symfun;
    end