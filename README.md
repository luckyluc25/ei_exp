# ei_exp

Those data and scripts underlie the publication : " Dissociable roles of cortical excitation-inhibition balance during patch-leaving versus value-guided decisions" by Kaiser, Gruendler, Luettgau, Speck and Jocham,2020.

If there is anything unclear in the code or if you have questions, please contact kaiserl[at]hhu.de

Anonymized data files can be found in the folder "logfiles". An anonymized table with all MRS results can be found in the folder 'ressources' as well as the extracted grey and white matter concentrations per participant. Due to ethical reasons, we cannot freely share the raw MR data.

All codes used to analyze data in the main manuscript are custom-written matlab codes. In order to analyze effect sizes you need to download the measures of effect size toolbox (https://github.com/hhentschke/measures-of-effect-size-toolbox) by Hentschke & Stüttgen (2018). All scripts used to analyze data in the patch-leaving phase can be found in the folder "First_Phase". All scripts used to analyze value-guided choice data can be found in the folder "Second_Phase". All remaining scripts are in the supplementary folder.

For the supplementary analysis, we used a combination of matlab and python 2 code (to fit the DDM models). Our python code is based on a repo published by Dr Anne Urai in 2019 (https://github.com/anne-urai/2019_Urai_choice-history-ddm; Urai, A. E., De Gee, J. W., Tsetsos, K., & Donner, T. H. (2019). Choice history biases subsequent evidence accumulation. Elife, 8, e46331.)
