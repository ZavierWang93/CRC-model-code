function funval = getOmega(vec_z, model_specs, contractor_idx)

for kkk=1:size(model_specs.omega{contractor_idx},1)
    omega = model_specs.omega{contractor_idx}(kkk,:);
    knot_x = log(1+model_specs.knot_omega{contractor_idx});

    knot_y = omega;

    funval(kkk,vec_z<=knot_x(1)) = knot_y(1);
    funval(kkk,vec_z>=knot_x(end)) = knot_y(end);
    fI = find(vec_z>knot_x(1) & vec_z<knot_x(end));
    funval(kkk,fI) = pchip(knot_x, knot_y, vec_z(fI));

    if (any(isnan(funval),'all'))
        error('Error in getOmega');
    end

end
