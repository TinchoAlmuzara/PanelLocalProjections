# PanelLocalProjections
This repository contains files to implement panel local projection estimation and inference.

Reference: [Almuzara Martin](https://martinalmuzara.com/research.html) and [Victor Sancibrian](https://sancibrian-v.github.io), 2024: “[Micro Responses to Macro Shocks.](https://www.newyorkfed.org/medialibrary/media/research/staff_reports/sr1090.pdf)” Federal Reserve Bank of New York Staff Report, no. 1090.

The <ins>**matlab_package**</ins> directory contains Matlab code:
  - The main function is **LP_panel.m**
  - It produces panel local projection estimates, standard errors, confidence intervals and p-values for zero response tests.
  - It can accommodate controls, different types of fixed effects, cumulative impules responses, etc.
  - It also implements the small-sample refinements suggested in the paper.
  - The file **usage_example.m** illustrates how to use it.

The <ins>**replication**</ins> directory contains files to replicate the simulation study of Section 4 and the empirical analysis of Section 5:
  - The empirical analysis uses Compustat and CRSP data which are proprietary and cannot be shared. This is why panel_data.csv appears empty in sec5_empirical_illustration/indata.  


