# PanelLocalProjections
This repository contains files to implement panel local projection estimation and inference.

Reference: [Almuzara Martin](https://martinalmuzara.com/research.html) and [Victor Sancibrian](https://sancibrian-v.github.io), 2024: “[Micro Responses to Macro Shocks.](https://www.newyorkfed.org/medialibrary/media/research/staff_reports/sr1090.pdf)” FRBNY Staff Report, no. 1090.

The <ins>**matlab_code**</ins> directory contains Matlab code:
  - The main function is **LP_panel.m**
  - It produces panel local projection estimates, standard errors, confidence intervals and p-values for zero response tests.
  - It can accommodate controls, different types of fixed effects, cumulative impules responses, etc.
  - It also implements the small-sample refinements suggested in the paper.
  - The file **usage_example.m** illustrates how to use it.

The <ins>**R_code**</ins> directory contains R code:
  - The main function is **LP_panel.R**
  - It produces panel local projection estimates, standard errors, confidence intervals and p-values for zero response tests.
  - It accommodates the same options and refinements as the Matlab code.
  - The file **usage_example.R** illustrates how to use it.


