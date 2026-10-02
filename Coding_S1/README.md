# Coding S1 — R analysis examples

These scripts are shared as adaptable examples for similar systematic-review and meta-analysis workflows. They are **not a package for reproducing the manuscript results**. Modify the input data, outcome definitions, comparison directions, pairing assumptions, dependence structure and model settings before using them for related primary analyses.

## Requirements and execution

Install R, then run in R:

    install.packages(c("metafor", "forestplot"))

From a terminal:

    Rscript Coding_S1/run_all.R

This runs the retained FoG example and saves plots, effect estimates, model summaries and session information under Coding_S1/outputs/. The source script contains example aggregate data. No separate extracted datasets or saved analysis outputs are distributed here.

The FoG example currently implements four instrument-specific analyses: questionnaire and assessment-course outcomes, each versus STN-only and baseline. It does not implement the manuscript overall FoG multilevel model. Only the output directory was changed from the designated source script; scripts/SMD_metaanlysis_source.R preserves its original contents.

For additional motor, medication, frequency and leave-one-study-out examples, first prepare your own CSV files in Coding_S1/data/ and run:

    Rscript Coding_S1/run_all.R --all

Required filenames: motor_overall.csv, motor_medication.csv, motor_high_risk_excluded.csv, frequency_extracted.csv. Modify scripts/analyze.R to fit your schema and analysis plan. It expects Research, Assessment, Medication, paired contrast means/SEs and alternative MD/SE fields for motor data, and frequency-specific final MD/SE fields for frequency data. Expected names are explicit in the code.

To implement similar **primary overall analyses**, modify the scripts to include the appropriate studies/outcomes and specify the correct model, sampling variances/covariances and study-level dependence. The provided stratified FoG example is not a substitute for an overall model. A within-subject correlation of 0.5 is used where assumed in the example; baseline correlations are recovered from usable change SDs when available. No sensitivity analysis varying correlation values is added.

## Files

- run_all.R: portable entry point.
- scripts/SMD_metaanlysis.R: runnable user-designated FoG example.
- scripts/SMD_metaanlysis_source.R: unchanged designated source.
- scripts/analyze.R: adaptable motor/subgroup/sensitivity examples; requires your own data.
- scripts/export_embedded_fog_data.R: exports aggregate FoG data embedded in the example.
- source_archive/: all retained historical R scripts. They may require editing, additional packages, input files or previously created R objects; they are not run by the entry point. Personal local paths are replaced with <ORIGINAL_LOCAL_PATH_REDACTED>.

Tested: portable FoG entry point with R 4.4.2, metafor 4.8-0 and forestplot. Historical fragments are supplied for reference and are not certified as independently runnable.

### Historical source file mapping

| Archived file | Original file |
|---|---|
| 01_original.R | 批量计算Mean difference.R |
| 02_original.R | 由中位数等计算均值.R |
| 03_original.R | Cross-over.R |
| 04_original.R | Frequency_meta.R |
| 05_original.R | Low and hIGH.R |
| 06_original.R | Mean difference的计算.R |
| 07_original.R | ON-OFF META.R |
| 08_original.R | SMD_metaanlysis.R |
| 09_original.R | forest plot\MD_difference_自动计算.R |
| 10_original.R | forest plot\ON_OFF_META_带prediciton.R |
| 11_original.R | forest plot\ON_OFF_META.R |
| 12_original.R | forest plot\sensitive_analysis\leave_one_out.R |
