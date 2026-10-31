function funval = logPriorSL(cohort,mc)

aa = 43.81;
bb = 155.56;
No_CS_pathway_pro=mc.No_CS_pathway_pro;

funval = sum((aa - 1) * log(cohort.age_specific_prop_CRC_SL*No_CS_pathway_pro) ...
    + (bb - 1) * log(1 - cohort.age_specific_prop_CRC_SL*No_CS_pathway_pro) - betaln(aa, bb));

end