ESPSS brings the New Statistics into SPSS (estimation, meta-analysis, sample-size planning, open science).

ESPSS is an extension bundle for SPSS.  It added custom menus and syntax to make it much easier to obtain the types of output required for the estimation approach (effect sizes and confidence intervals).  It *de-emphasizes* the null-hypothesis significance testing (NHST) approach.

Commands currently supported:
* ETTEST - estimates the difference between two independent means (an independent samples t-test).  Compared to the standard t-test in SPSS, this command provides Cohen's d and its confidence interval, and the CIs for each group mean.  It also emphasizes the effect size (Mdiff) and its CI rather than t and p.
* EPAIRED - estimates the difference between two related means (a paired-samples t-test aka a repeated-measures t-test).  Compared to the standard paired t-test in SPSS, this command provides COhen's d and ints confidence interval and the CIs for each group mean.  It also emphasizes the effect size (Mdiff) and its CI rather than t and p.

ESPSS was written by Bob Calin-Jageman.  Please post questions/bug reports to the GitHub page:
* https://github.com/rcalinjageman/ESPSS

Roadmap:
* First, complete estimation-approach versions of the basic stats toolkit:
** Simple linear correlation
** One-sample t-test
** Descriptives
** Frequencies
** Comparison of two independent proportions (what is Chi square in NHST approach)
** One-way ANOVA
** Repeated-measures and/or factorial ANOVA
* Then, see if it is possible to add some reasonable visualizations

Principles for ESPSS are:
* Try to make things easy on the user
* Try to stick as close to SPSS in terms of syntax and UI as possible (though SPSS doesn't make this easy)
* Re-use existing SPSS commands/output wherever possible rather than re-inventing the wheel


ESPSS Contains Code Lifted from Others.  I've tried to obtain correct permissions wherever possible:
* NONCT2.sps by M. J. Smithson's site. with edits by K. Wuensch
** http://core.ecu.edu/psyc/wuenschk/SPSS/SPSS-Programs.htm
* The SPSSwrap function from spssaux2.py by JKP, IBM
** https://www.ibm.com/developerworks/community/files/app#/file/108ce263-b14e-4375-9f70-6d1f733d132b
* STATS_CORRELATION.py by ___
**
* Code adapted from ESCI by Geoff Cumming
* And more