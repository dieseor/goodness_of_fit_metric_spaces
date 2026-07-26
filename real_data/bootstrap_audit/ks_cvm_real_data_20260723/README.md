# KS/CvM bootstrap audit for selected real-data cases

The table reports the conditional Monte Carlo uncertainty of each bootstrap p-value. `bootstrap_sd` is the standard deviation of the statistic over bootstrap replicates; `mc_se` is the Monte Carlo standard error of the tail probability estimate, not a sampling standard error for the data-generating model.

## Bootstrap-statistic summaries

```text
                           case          method statistic    B   observed
                  pogojump_fast fast_multiplier        KS 5000 1.13389342
                  pogojump_fast fast_multiplier       CVM 5000 0.06832806
           pogojump_reestimated     reestimated        KS 5000 1.13389342
           pogojump_reestimated     reestimated       CVM 5000 0.06832806
 risoe_nov_dec_125m_start4_fast fast_multiplier        KS 5000 1.07787533
 risoe_nov_dec_125m_start4_fast fast_multiplier       CVM 5000 0.06515016
 exceedances   p_value       mc_se mc_wilson_lower mc_wilson_upper
         489 0.0979804 0.004200837      0.08987265       0.1063449
        2614 0.5228954 0.007063712      0.50894318       0.5366218
         628 0.1257748 0.004686676      0.11670074       0.1350741
        2571 0.5142971 0.007068216      0.50034096       0.5280372
        2156 0.4313137 0.007003807      0.41753078       0.4449749
        3676 0.7352529 0.006239887      0.72279284       0.7472460
 bootstrap_mean bootstrap_sd  observed_z lag1_correlation
     0.87358170   0.19459092  1.33773827     -0.021116082
     0.07951881   0.03895050 -0.28730705     -0.010007300
     0.90375466   0.19881845  1.15753220     -0.011575418
     0.07771688   0.03639332 -0.25798213      0.003718001
     1.06616299   0.18043734  0.06491086      0.015898484
     0.08675904   0.03192617 -0.67683936      0.020112534
```

## Stored fast versus stored fully re-estimated multiplier bootstrap (PogoJump)

```text
 statistic p_value_fast  mc_se_fast observed_fast bootstrap_mean_fast
       CVM    0.5228954 0.007063712    0.06832806          0.07951881
        KS    0.0979804 0.004200837    1.13389342          0.87358170
 bootstrap_sd_fast p_value_reestimated mc_se_reestimated observed_reestimated
         0.0389505           0.5142971       0.007068216           0.06832806
         0.1945909           0.1257748       0.004686676           1.13389342
 bootstrap_mean_reestimated bootstrap_sd_reestimated p_value_difference
                 0.07771688               0.03639332         0.00859828
                 0.90375466               0.19881845        -0.02779444
 combined_mc_se difference_in_mc_se_units
    0.009992783                 0.8604491
    0.006293803                -4.4161598
```

Interpretation notes:

- A large difference between the calibrated KS and CvM p-values is not, by itself, an error: KS is a supremum functional, while CvM integrates squared departures. Their null distributions are different, so their p-values admit no general ordering.
- The `p_reconstruction_error` column must be zero up to floating-point precision whenever the stored raw replicates are available.
- The Risoe run is saved in this directory because the original paper-table CSV retained only p-values, not the bootstrap-statistic vectors.
