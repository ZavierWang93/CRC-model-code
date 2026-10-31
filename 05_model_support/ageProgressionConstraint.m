function fun_val=ageProgressionConstraint(model_specs)

fun_val = 0;
age_effect = model_specs.age_effect_on_disease_progression;

diff_v = diff([1 age_effect]);

if any(diff_v < 0)
    fun_val=-100000;
end

end