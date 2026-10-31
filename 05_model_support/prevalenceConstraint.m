function fun_val=prevalenceConstraint(model_specs,mc)

fun_val = 0;
prevalence = model_specs.prevalence;
Prevalence_constraint=mc.Prevalence_constraint;

if (Prevalence_constraint == 1)

    isSmaller = all(prevalence(end) < prevalence(1:(end-1)));
    if   ~isSmaller
        fun_val = -100000;
    end

end

end