function funval = logPriorNA(ms_final)

    alpha_NA = 2548;
    beta_NA = 61153;
    NA_p = ms_final.exit_rate(:,1);

    log_prior_NA= (alpha_NA - 1) * log(NA_p) ...
                        + (beta_NA - 1) * log(1 - NA_p);

    funval = log_prior_NA;

end