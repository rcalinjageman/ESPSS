# -*- coding: utf-8 -*-

"""ETTEST extension command"""

__author__ =  'RC-J'
__version__=  '1.0.0'

import spss, spssaux, spssdata
from extension import Template, Syntax, processcmd
import random
from xml.etree import ElementTree
import re
import math
import textwrap

# This it the python module for an extension that completes a t-test with an estimation approach (ETTEST)
# The improvements are: reporting cohen's d and the CI for cohen's d.
# Cohen's d is reported correct for bias (dunbiased = Hedges g)

# This module was written by Bob Calin-Jageman

# This module contains code from:
#       NONCT2.sps by M. J. Smithson's site.
#       Edits by K. Wuensch
#           http://core.ecu.edu/psyc/wuenschk/SPSS/SPSS-Programs.htm
#       The SPSSwrap function from spssaux2.py by JKP, IBM
#           https://www.ibm.com/developerworks/community/files/app#/file/108ce263-b14e-4375-9f70-6d1f733d132b
#       And lots of inspiration from the extension STATS_CORRELATION.py

# To do:
#   Ensure correct function with different error conditions (no data, no levels, etc.)
#   Make options portion of syntax optional (not sure why it is not)
#   Add cut-point option
#   Add graphical output
#   Find a more elegant solution that just re-using the syntax for Cohen's d CI; re-implement in python
#   Find a better way to seperate the syntax, output text, python logic to make it easier to maintain this code
#   Find out how to integrate current version into output note


def do_ettest(groups, gone, gtwo, dvlist, missing, conf):
    #This is the actual function that completes the ETTEST procedure
      #It accepts the list of arguments parsed out in the run subroutine SPSS initially calls
      #The run subroutine calls this function via SPSS' processcmd routine

    #Set to True to get debugging output; set to False for distribution
    debugit = False
    
    if debugit:
        print("Entered do_ettest")
        print(groups)
        print(gone)
        print(gtwo)
        print(dvlist)
        print missing
        print conf

    #Get the name for the current dataset
    activeds = spss.ActiveDataset()

    #Pre-req checks - Make sure all requirements are met to run the command
    #To do: write this section


    #Put the parameters for the syntax command back into string format
    groups_string = " ".join(groups)
    dvlist_string = " ".join(dvlist)
    gone_string = str(gone)
    gtwo_string = str(gtwo)
  
    #Setup the t-test command--will need to adapt this for the more complex t-test commands possible
    mycommand="""T-TEST GROUPS=%(groups_string)s(%(gone_string)s %(gtwo_string)s)
/MISSING=ANALYSIS
/VARIABLES=%(dvlist_string)s
/CRITERIA=CI(0.%(conf)s).""" % locals()  
    if debugit:
        print(mycommand)
       
    #Create the descriptive stats output
    tdestable,err = spssaux.CreateDatasetOutput(mycommand, omsid='T-Test', subtype='Group Statistics')
    #Probably should do some checking on the spss failcode as well    
    if not err == 0:
        raise ValueError(_("Descriptive Stats Command Failed"))
    spss.Submit("DATASET ACTIVATE "+tdestable+".")

    #Now collect the data in an array of t-test results
    desdata = spssdata.Spssdata()
    resultlist = []
    start = True
    count = 0

    for row in desdata:
        if start:
            resultlist.append(tresult())
            resultlist[count].conf = conf
            resultlist[count].iv = groups
            resultlist[count].dv = str(row.Var1)
            resultlist[count].level1 = str(row.Var2)
            resultlist[count].level1_code = gone
            resultlist[count].level2_code = gtwo
            resultlist[count].n1 = row.N
            resultlist[count].m1 = row.Mean
            resultlist[count].s1 = row[7]
            resultlist[count].sem1 = row[8]
            #Would be nice if we 
            start = False
        else:
            resultlist[count].conf = conf
            resultlist[count].level2 = str(row.Var2)
            resultlist[count].n2 = row.N
            resultlist[count].m2 = row.Mean
            resultlist[count].s2 = row[7]
            resultlist[count].sem2 = row[8]
            start = True
            count = count + 1
            
    desdata.CClose()
    spss.Submit("DATASET ACTIVATE "+activeds+".")
    spss.Submit("DATASET CLOSE "+tdestable+".")

    #Now do the t-table output
    ttable,err = spssaux.CreateDatasetOutput(mycommand, omsid='T-Test', subtype='Independent Samples Test')
    #Probably should do some checking on the spss failcode as well    
    if not err == 0:
        raise ValueError(_("T-TEST Command Failed"))
    spss.Submit("DATASET ACTIVATE "+ttable+".")
    
    #Now update the t-table output to have the columns needed for the CI of d script
    mydcisyntax = dcisyntax()
    spss.Submit(mydcisyntax.create_data_syntax)
    spss.Submit("COMPUTE conf=0."+str(conf)+""".
EXECUTE.""")

    #Now gather the t-table output while also writing the sample sizes into the output file
    tdata = spssdata.Spssdata(accessType='w')
    tdata.append(spssdata.vdef('n1', vfmt=["F",8,0]))
    tdata.append(spssdata.vdef('n2', vfmt=["F",8,0]))
    tdata.commitdict()
    countt = 0
    start = True

    for row in tdata:
        if start:
            resultlist[countt].eq_t = row.t
            resultlist[countt].eq_df = row.df
            resultlist[countt].eq_p = row[9]
            resultlist[countt].mdiff = row.MeanDifference
            resultlist[countt].eq_mdiff_sed = row[11]
            resultlist[countt].eq_mdiff_low = row.Lower
            resultlist[countt].eq_mdiff_high = row.Upper
            tdata.casevalues([resultlist[countt].n1, resultlist[countt].n2])
            start = False
        else:
            resultlist[countt].ne_t = row.t
            resultlist[countt].ne_df = row.df
            resultlist[countt].ne_p = row[9]
            resultlist[countt].ne_mdiff_sed = row[11]
            resultlist[countt].ne_mdiff_low = row.Lower
            resultlist[countt].ne_mdiff_high = row.Upper
            tdata.casevalues([resultlist[countt].n1, resultlist[countt].n2])
            start = True
            countt = countt + 1

    tdata.CClose()
    

    #we are ready to use the pre-existing ci for d script
    if debugit:
        print "Ready to run ci on d syntax"


    tw = SpssWrap()
    cmdlist = tw.wrap(mydcisyntax.do_d_ci_syntax)
    spss.Submit(cmdlist)
    
    
    if debugit:
        print "Completed ci on d syntax"

    #We need to calculate critical t for making individual mean CIs and difference MoE
    
    ci_level = float(conf)/100
    ci_level = ci_level + ((1-ci_level)/2)
    
    gettcrit = """COMPUTE m1tcrit = IDF.T("""+str(ci_level)+""", n1-1).
COMPUTE m2tcrit =IDF.T("""+str(ci_level)+""", n2-1).
COMPUTE mdifftcrit =IDF.T("""+str(ci_level)+""", df).
EXECUTE."""
    spss.Submit(gettcrit)


    #one more cycle through the data to collect the cis for d
    tdata = spssdata.Spssdata()
    countt = 0
    start = True

    for row in tdata:
        if start:
            resultlist[countt].eq_d_low  = row.lowd
            resultlist[countt].eq_d_high = row.highd
            resultlist[countt].m1_tcrit = row.m1tcrit
            resultlist[countt].m2_tcrit = row.m2tcrit
            resultlist[countt].eq_mdiff_tcrit = row.mdifftcrit
            if resultlist[countt].eq_t < 0:
                #If t was negative, need to invert CI
                tempvar = resultlist[countt].eq_d_low * -1
                resultlist[countt].eq_d_low = resultlist[countt].eq_d_high * -1
                resultlist[countt].eq_d_high = tempvar
            start = False
        else:
            resultlist[countt].ne_d_low  = row.lowd
            resultlist[countt].ne_d_high = row.highd
            resultlist[countt].ne_mdiff_tcrit = row.mdifftcrit
            if resultlist[countt].ne_t < 0:
                #If t was negative, need to invert CI
                tempvar = resultlist[countt].ne_d_low * -1
                resultlist[countt].ne_d_low = resultlist[countt].ne_d_high * -1
                resultlist[countt].ne_d_high =  tempvar
            start = True
            if debugit:
                resultlist[countt].printout()
            countt = countt + 1
    tdata.CClose()
    
    spss.Submit("DATASET ACTIVATE "+activeds+".")
    spss.Submit("DATASET CLOSE "+ttable+".")

    #Don't forget to patch up the results from SPSS, filling in missing info (CI on d, MoE and CI for individual means, MoE for MDiff
    for myresult in resultlist:
            myresult.calc_missing()
  
    #Ok - we are ready for reporting!
    display(resultlist)
 
    if debugit:
        print("Current end of do_ettest")
    #That's the end of ettest procedure

def display(resultlist):
    #Do the actual output of the e t-test
    spss.StartProcedure(_("Estimation: Compare two groups"), "ETTEST")

    for myresult in resultlist:

        myrowlabels = [myresult.level1, myresult.level2, "Raw Score Difference", "Standardized difference (Cohen's d)"]
        mycolumnlabels = ["M", str(myresult.conf)+"% MoE", str(myresult.conf)+"% CI [low, high]", "N", "s"]
        m1ci = "[{0:.3f}, {1:.3f}]".format(myresult.m1_low, myresult.m1_high)
        m2ci = "[{0:.3f}, {1:.3f}]".format(myresult.m2_low, myresult.m2_high)

        if myresult.varequal():
            dci = "[{0:.3f}, {1:.3f}]".format(myresult.eq_d_low, myresult.eq_d_high)
            mci = "[{0:.3f}, {1:.3f}]".format(myresult.eq_mdiff_low, myresult.eq_mdiff_high)
            mycells = [[myresult.m1, myresult.m1_moe, m1ci, myresult.n1, myresult.s1],
                     [myresult.m2, myresult.m2_moe, m2ci, myresult.n2, myresult.s2],
                     [myresult.mdiff, myresult.eq_mdiff_moe, mci, "-", "-"],
                     [myresult.d, "-", dci, "-", "-"]
                     ]
            mycaption = "Equal variance assumed because group standard deviations are similar (within a factor of 2 of each other)."

            mytable = spss.BasePivotTable("Comparison for " + myresult.dv, "ETTEST", caption=mycaption)
            mytable.SimplePivotTable(rowdim = "Row",
                                      rowlabels= myrowlabels,
                                      coldim = "Column",
                                      collabels = mycolumnlabels,
                                      cells = mycells
                                      )
            eq_tresult =  "NHST Approach for " + myresult.dv + ": t({0:.0f}) = {1:.2f}, p = {2:.4f}, without equal variance assumed.".format(myresult.eq_df, myresult.eq_t, myresult.eq_p)
            eq_txtblock = spss.TextBlock("NHST Approach for "+ myresult.dv, eq_tresult)
        else:
            nedci = "[{0:.3f}, {1:.3f}]".format(myresult.ne_d_low, myresult.ne_d_high)
            nemci = "[{0:.3f}, {1:.3f}]".format(myresult.ne_mdiff_low, myresult.ne_mdiff_high)
            necells = [[myresult.m1, myresult.m1_moe, m1ci, myresult.n1, myresult.s1],
                     [myresult.m2, myresult.m2_moe, m2ci, myresult.n2, myresult.s2],
                     [myresult.mdiff, myresult.ne_mdiff_moe, nemci, "-", "-"],
                     [myresult.d, "-", nedci, "-", "-"]
                     ]
            mycaption = "Results reported without assuming equal variance because group standard deviations are not similar (greater than 2-fold difference).  CI for d is not calculated in this case (working on it)."

            netable = spss.BasePivotTable("Comparison for " + myresult.dv, "ETTEST", caption=mycaption)
            netable.SimplePivotTable(rowdim = "Row",
                                      rowlabels= myrowlabels,
                                      coldim = "Column",
                                      collabels = mycolumnlabels,
                                      cells = necells
                                      ) 
            ne_tresult =  "NHST Approach for " + myresult.dv + ": t({0:.2f}) = {1:.2f}, p = {2:.4f}, without the assumption of equal variance.".format(myresult.ne_df, myresult.ne_t, myresult.ne_p)
            ne_txtblock = spss.TextBlock("NHST Approach "+ myresult.dv, ne_tresult)

    othercomments = spss.TextBlock("Important Notes", """

The ETTEST routine was written by Bob Calin-Jageman, v1.0, 4/26/2017
You will see a series of "Warning #534 messages in the notes for this command.
These messages are normal and can be ignored.

Some tips:
* Visualize your data--don't draw conclusions until you've made a good visualization showing
all the data that emphasizes the group comparison.  Try ESCI for this.
* Don't draw causal conclusions unless random assignment was used.
* The confidence intervals only make sense to the degree to which the sample
is representative of the population.
* Clearly predicted (confirmatory) results are to be trusted more than those found by
exploring the data, but in either case seek replication.
* Be sure to present the whole story when reporting results; statistics can only be judged in
the context of all the data collected and analyses run.
""")
        
    spss.EndProcedure()
    
class tresult:
    #This is a class to hold the results of a t-test
    # It's not very elegant, partly because I'm not a python whiz
    # Partly because it's a messy structure.  Results for both equal variance assumed and not assumed
    #   Are mashed together into this same structure..just felt a little easier, though not as abstract
    #   Could use some improvement, for sure
    conf = 0
    n1 = 0
    n2 = 0
    m1 = 0
    m2 = 0
    m1_tcrit = 0
    m2_tcrit = 0
    m1_low = 0
    m1_high = 0
    m2_low = 0
    m2_high = 0
    m1_moe = 0
    m2_moe = 0
    s1 = 0
    s2 = 0
    sem1 = 0
    sem2 = 0
    d = 0
    eq_d_low = 0
    eq_d_high = 0
    ne_d_low = 0
    ne_d_high = 0
    mdiff = 0
    eq_mdiff_moe = 0
    eq_mdiff_sed = 0
    eq_mdiff_low = 0
    eq_mdiff_high = 0
    eq_mdiff_tcrit = 0
    ne_mdiff_moe = 0
    ne_mdiff_sed = 0
    ne_mdiff_low = 0
    ne_mdiff_high = 0
    ne_mdiff_tcrit = 0
    iv = "IV"
    dv = "DV"
    level1 = "Level1"
    level2 = "Level2"
    level1_code = 0
    level2_code = 0
    eq_t = 0
    eq_df = 0
    eq_p = 0
    ne_t = 0
    ne_df = 0
    ne_p = 0
    tbl = 0

    def printout(self):
        print "Compare two groups report"
        print "IV = ", self.iv, " DV = ", self.dv
        for property, value in vars(self).iteritems():
            print property, ": ", value
        print "assume equal variance: d = ", self.d, "95% CI [", self.eq_d_low, ", ", self.eq_d_high, "]"
        print "no assume equal variance: d = ", self.d, "95% CI [", self.ne_d_low, ", ", self.ne_d_high, "]"
        print "End report"

    def calc_missing(self):
        #Calculate cohen's d from t, n1, and n2.
        #Adjusted for small upward bias using formula reported by Lakens, 2013
        self.d = self.eq_t * math.sqrt(self.n1+self.n2)/math.sqrt(self.n1*self.n2)
        self.d = self.d * (1 - (3 / (4 * (self.n1+self.n2) - 9)) )

        #Calculate mean difference and the raw score CI
        self.mdiff = self.m2 - self.m1
        self.eq_mdiff_moe = self.eq_mdiff_tcrit*self.eq_mdiff_sed
        #self.eq_mdiff_low = self.mdiff - self.eq_mdiff_moe
        #self.eq_mdiff_high = self.mdiff + self.eq_mdiff_moe
        self.ne_mdiff_moe = self.ne_mdiff_tcrit*self.ne_mdiff_sed
        #self.ne_mdiff_low = self.mdiff - self.ne_mdiff_moe
        #self.ne_mdiff_high = self.mdiff + self.ne_mdiff_moe

        #Calculate individual mean CIs
        self.m1_moe = self.m1_tcrit * self.sem1
        self.m1_low = self.m1 - self.m1_moe
        self.m1_high = self.m1 + self.m1_moe
        self.m2_moe = self.m2_tcrit * self.sem2
        self.m2_low = self.m2 - self.m2_moe
        self.m2_high = self.m2 + self.m2_moe

        #clean up some strings
        self.dv = str(self.dv).strip()
        self.iv = str(self.iv).strip()
        self.level1 = str(self.level1).strip()
        self.level2 = str(self.level2).strip()

    def varequal(self):
        if self.s1 > (self.s2 * 2):
            return False
        elif self.s2 > (self.s1 * 2):
            return False
        else:
            return True


class dcisyntax:
    #This is a class to hold the extensive syntax required to caluclate the CI on cohen's d
    #The syntax is from MJ Smithson with comments and tewaks from K. Wuensch
    #Wuensch warns that this syntax doesn't always handle negative t values properly
    #So this module always uses abs(t) and then, once the CI is calculated, inverts the CI if t<0
    #The initial create_data_syntax string creates the additional columns needed for the script
    #If you want a different confidence level, need to re-calculate conf after calling create_data_syntax
    #The main syntax is in do_d_ci_syntax.  It is too long to be added directly with submit,
    #Instead it needs to be submitted using SpssWrap, below
    create_data_syntax = """
COMPUTE tval=abs(t).
COMPUTE conf=0.95.
NUMERIC lc2 ucdf uc2 lcdf power d lowd highd (F8.4).
EXECUTE.
"""

    do_d_ci_syntax = """
COMMENT	From NONCT2.sps at M. J. Smithson's site.

COMMENT	THIS SCRIPT COMPUTES CONFIDENCE INTERVALS FOR THE NONCENTRALITY PARAMETER
COMMENT	FOR THE NONCENTRAL T DISTRIBUTION.
COMMENT	IT USES THE SPSS NONCENTRAL T CALCULATOR AND LAUBSCHER'S (1960) NORMAL APPROXIMATION
COMMENT	TO THE NONCENTRAL F WITH 1 DF FOR THE SPECIAL CASE WHERE F = T^2,
COMMENT	WITH A DECISION RULE FOR CHOOSING BETWEEN THEM. THE REASON FOR THIS IS THAT THE 
COMMENT	NONCENTRAL T ALGORITHM FAILS FOR LARGE SAMPLE SIZE OR EFFECT SIZE.
COMMENT	THE FIRST PART USES THE NONCENTRAL T CALCULATOR IN SPSS.
COMMENT	THIS COMPUTES THE LOWER LIMIT ON THE T STATISTIC.
COMPUTE #LC3 = TVAL .
COMPUTE LC2 = TVAL/2 .
COMPUTE #LC1 = -TVAL .
COMPUTE #CUMF1 = NCDF.T(TVAL,DF,#LC1) .
COMPUTE #ULIM = 1-(1-CONF)/2 .
LOOP IF (#CUMF1 LT #ULIM) .
+	COMPUTE LC2 = #LC1 .
+	COMPUTE #LC1 = #LC1 - TVAL .
+	COMPUTE #CUMF1 = NCDF.T(TVAL,DF,#LC1) .
END LOOP .
COMPUTE #CUMF2 = NCDF.T(TVAL,DF,LC2) .
COMPUTE #DIFF = 1 .
LOOP IF (#DIFF GT .00005) .
+	DO IF (#CUMF2 LT #ULIM) .
+		COMPUTE #LC3 = LC2 .
+		COMPUTE LC2 = (LC2 + #LC1)/2 .
+		COMPUTE #CUMF2 = NCDF.T(TVAL,DF,LC2) .
+	ELSE .
+		COMPUTE #LC1 = LC2 .
+		COMPUTE LC2 = (LC2 + #LC3)/2 .
+		COMPUTE #CUMF2 = NCDF.T(TVAL,DF,LC2) .
+	END IF .
+	COMPUTE #DIFF = ABS(#CUMF2 - #ULIM) .
END LOOP .
COMPUTE UCDF  = NCDF.T(TVAL,DF,LC2) .
EXECUTE .
COMMENT	
COMMENT	THIS COMPUTES THE UPPER LIMIT ON THE T STATISTIC.
COMPUTE #UC3 = 2*TVAL .
COMPUTE UC2 = 1.5*TVAL .
COMPUTE #UC1 = TVAL .
COMPUTE #LLIM = (1-CONF)/2 .
COMPUTE #CUMF3 = NCDF.T(TVAL,DF,#UC3) .
LOOP IF (#CUMF3 GT #LLIM) .
+	COMPUTE UC2 = #UC3 .
+	COMPUTE #UC3 = #UC3 + TVAL .
+	COMPUTE #CUMF3 = NCDF.T(TVAL,DF,#UC3) .
END LOOP .
COMPUTE #CUMF2 = NCDF.T(TVAL,DF,UC2) .
COMPUTE #DIFF = 1 .
LOOP IF (#DIFF GT .00001) .
+	DO IF (#CUMF2 LT #LLIM) .
+		COMPUTE #UC3 = UC2 .
+		COMPUTE UC2 = (UC2 + #UC1)/2 .
+		COMPUTE #CUMF2 = NCDF.T(TVAL,DF,UC2) .
+	ELSE .
+		COMPUTE #UC1 = UC2 .
+		COMPUTE UC2 = (UC2 + #UC3)/2 .
+		COMPUTE #CUMF2 = NCDF.T(TVAL,DF,UC2) .
+	END IF .
+	COMPUTE #DIFF = ABS(#CUMF2 - #LLIM) .
END LOOP .
COMPUTE LCDF  = NCDF.T(TVAL,DF,UC2) .
COMMENT	
COMMENT	THIS NEXT STATEMENT COMPUTES THE POWER IN RELATION TO THE T VALUE.
COMPUTE POWER = 1 - NCDF.T(IDF.T(1-(1-CONF)/2,DF),DF,TVAL) .
EXECUTE .
COMMENT	
COMMENT	THE SECOND PART USES LAUBSCHER'S SQUARE-ROOT APPROXIMATION.
COMMENT	THIS COMPUTES THE LOWER LIMIT ON THE F NONCENTRALITY PARAMETER.
COMPUTE #LLC3 = TVAL**2 .
COMPUTE #LLC2 = TVAL**2/2 .
COMPUTE #LLC1 = .001 .
COMPUTE #ULIM = 1-(1-CONF)/2 .
COMPUTE #CUMF1 = 1-CDFNORM((Sqrt(2*(1+#LLC1)-((1+2*#LLC1)/(1+#LLC1)))-Sqrt((2*DF-1)*TVAL**2*1/DF))/
  Sqrt(1*TVAL**2/DF+((1+2*#LLC1)/(1+#LLC1)))) .
LOOP IF (#CUMF1 LT #ULIM) .
+	COMPUTE #LLC2 = #LLC1 .
+	COMPUTE #LLC1 = #LLC1/4 .
+	COMPUTE #CUMF1 = 1-CDFNORM((Sqrt(2*(1+#LLC1)-((1+2*#LLC1)/(1+#LLC1)))-Sqrt((2*DF-1)*TVAL**2*1/DF))/
  Sqrt(1*TVAL**2/DF+((1+2*#LLC1)/(1+#LLC1)))) .
END LOOP .
COMPUTE #CUMF3 = 1-CDFNORM((Sqrt(2*(1+#LLC3)-((1+2*#LLC3)/(1+#LLC3)))-Sqrt((2*DF-1)*TVAL**2*1/DF))/
  Sqrt(1*TVAL**2/DF+((1+2*#LLC3)/(1+#LLC3)))) .
LOOP IF (#CUMF3 GT #ULIM) .
+	COMPUTE #LLC2 = #LLC3 .
+	COMPUTE #LLC3 = #LLC3 + TVAL**2 .
+	COMPUTE #CUMF3 = 1-CDFNORM((Sqrt(2*(1+#LLC3)-((1+2*#LLC3)/(1+#LLC3)))-Sqrt((2*DF-1)*TVAL**2*1/DF))/
  Sqrt(1*TVAL**2/DF+((1+2*#LLC3)/(1+#LLC3)))) .
END LOOP .
COMPUTE #CUMF2 = 1-CDFNORM((Sqrt(2*(1+#LLC2)-((1+2*#LLC2)/(1+#LLC2)))-Sqrt((2*DF-1)*TVAL**2*1/DF))/
  Sqrt(1*TVAL**2/DF+((1+2*#LLC2)/(1+#LLC2)))) .
COMPUTE #DIFF = 1 .
LOOP IF (#DIFF GT .00005) .
+	DO IF (#CUMF2 LT #ULIM) .
+		COMPUTE #LLC3 = #LLC2 .
+		COMPUTE #LLC2 = (#LLC2 + #LLC1)/2 .
+		COMPUTE #CUMF2 = 1-CDFNORM((Sqrt(2*(1+#LLC2)-((1+2*#LLC2)/(1+#LLC2)))-Sqrt((2*DF-1)*TVAL**2*1/DF))/
  Sqrt(1*TVAL**2/DF+((1+2*#LLC2)/(1+#LLC2)))) .
+	ELSE .
+		COMPUTE #LLC1 = #LLC2 .
+		COMPUTE #LLC2 = (#LLC2 + #LLC3)/2 .
+		COMPUTE #CUMF2 = 1-CDFNORM((Sqrt(2*(1+#LLC2)-((1+2*#LLC2)/(1+#LLC2)))-Sqrt((2*DF-1)*TVAL**2*1/DF))/
  Sqrt(1*TVAL**2/DF+((1+2*#LLC2)/(1+#LLC2)))) .
+	END IF .
+	COMPUTE #DIFF = ABS(#CUMF2 - #ULIM) .
END LOOP .
COMPUTE #UUCDF  = 1-CDFNORM((Sqrt(2*(1+#LLC2)-((1+2*#LLC2)/(1+#LLC2)))-Sqrt((2*DF-1)*TVAL**2*1/DF))/
  Sqrt(1*TVAL**2/DF+((1+2*#LLC2)/(1+#LLC2)))) .
COMMENT	
COMMENT	THIS COMPUTES THE UPPER LIMIT ON THE T NONCENTRALITY PARAMETER.
COMPUTE #UUC3 = 3*TVAL**2 .
COMPUTE #UUC2 = 2*TVAL**2 .
COMPUTE #UUC1 = TVAL**2 .
COMPUTE #LLIM = (1-CONF)/2 .
COMPUTE #CUMF1 = 1-CDFNORM((Sqrt(2*(1+#UUC1)-((1+2*#UUC1)/(1+#UUC1)))-Sqrt((2*DF-1)*TVAL**2*1/DF))/
  Sqrt(1*TVAL**2/DF+((1+2*#UUC1)/(1+#UUC1)))) .
LOOP IF (#CUMF1 LT #LLIM) .
+	COMPUTE #UUC2 = #UUC1 .
+	COMPUTE #UUC1 = #UUC1/4 .
+	COMPUTE #CUMF1 = 1-CDFNORM((Sqrt(2*(1+#UUC1)-((1+2*#UUC1)/(1+#UUC1)))-Sqrt((2*DF-1)*TVAL**2*1/DF))/
  Sqrt(1*TVAL**2/DF+((1+2*#UUC1)/(1+#UUC1)))) .
END LOOP .
COMPUTE #CUMF3 = 1-CDFNORM((Sqrt(2*(1+#UUC3)-((1+2*#UUC3)/(1+#UUC3)))-Sqrt((2*DF-1)*TVAL**2*1/DF))/
  Sqrt(1*TVAL**2/DF+((1+2*#UUC3)/(1+#UUC3)))) .
LOOP IF (#CUMF3 GT #LLIM) .
+	COMPUTE #UUC2 = #UUC3 .
+	COMPUTE #UUC3 = #UUC3 + TVAL**2 .
+	COMPUTE #CUMF3 = 1-CDFNORM((Sqrt(2*(1+#UUC3)-((1+2*#UUC3)/(1+#UUC3)))-Sqrt((2*DF-1)*TVAL**2*1/DF))/
  Sqrt(1*TVAL**2/DF+((1+2*#UUC3)/(1+#UUC3)))) .
END LOOP .
COMPUTE #CUMF2 = 1-CDFNORM((Sqrt(2*(1+#UUC2)-((1+2*#UUC2)/(1+#UUC2)))-Sqrt((2*DF-1)*TVAL**2*1/DF))/
  Sqrt(1*TVAL**2/DF+((1+2*#UUC2)/(1+#UUC2)))) .
COMPUTE #DIFF = 1 .
LOOP IF (#DIFF GT .00001) .
+	DO IF (#CUMF2 LT #LLIM) .
+		COMPUTE #UUC3 = #UUC2 .
+		COMPUTE #UUC2 = (#UUC2 + #UUC1)/2 .
+		COMPUTE #CUMF2 = 1-CDFNORM((Sqrt(2*(1+#UUC2)-((1+2*#UUC2)/(1+#UUC2)))-Sqrt((2*DF-1)*TVAL**2*1/DF))/
  Sqrt(1*TVAL**2/DF+((1+2*#UUC2)/(1+#UUC2)))) .
+	ELSE .
+		COMPUTE #UUC1 = #UUC2 .
+		COMPUTE #UUC2 = (#UUC2 + #UUC3)/2 .
+		COMPUTE #CUMF2 = 1-CDFNORM((Sqrt(2*(1+#UUC2)-((1+2*#UUC2)/(1+#UUC2)))-Sqrt((2*DF-1)*TVAL**2*1/DF))/
  Sqrt(1*TVAL**2/DF+((1+2*#UUC2)/(1+#UUC2)))) .
+	END IF .
+	COMPUTE #DIFF = ABS(#CUMF2 - #LLIM) .
END LOOP .
COMPUTE #LLCDF  = 1-CDFNORM((Sqrt(2*(1+#UUC2)-((1+2*#UUC2)/(1+#UUC2)))-Sqrt((2*DF-1)*TVAL**2*1/DF))/
  Sqrt(1*TVAL**2/DF+((1+2*#UUC2)/(1+#UUC2)))) .
COMMENT	
COMMENT	THIS NEXT STATEMENT COMPUTES THE POWER IN RELATION TO THE F VALUE.
COMPUTE #PPOWER = CDFNORM((Sqrt(2*(1+TVAL**2*(1/DF)*(1+DF+1))-((1+2*TVAL**2*(1/DF)*(1+DF+1))/(1+TVAL**2*(1/DF)*(1+DF+1))))-Sqrt((2*DF-1)*IDF.F(1-(1-CONF)/2,1,DF)*1/DF))/
  Sqrt(1*IDF.F(1-(1-CONF)/2,1,DF)/DF+((1+2*TVAL**2*(1/DF)*(1+DF+1))/(1+TVAL**2*(1/DF)*(1+DF+1))))) .
COMMENT	
COMMENT	NOW CHOOSE THE METHOD TO USE FOR THE FINAL ESTIMATES. 
COMMENT	THIS DECISION IS BASED ON THE OBSERVATION THAT SPSS NONCENTRAL T 
COMMENT	PRODUCES UPPER BOUND ESTIMATES THAT FALL BELOW ACCURATE VALUES UNTIL
COMMENT	ARE SIMPLY EQUAL TO THE TVAL THAT IS GIVEN AS THE STARTING POINT. 
COMMENT	THE NORMAL APPROXIMATION, ON THE OTHER HAND, TENDS TO UNDERESTIMATE 
COMMENT	THE TRUE NONCENTRALITY PARAMETER FOR SMALL N AND/OR SMALL EFFECTS. 
COMMENT	NONCENTRAL T ALSO PRODUCES A LOWER BOUND THAT HITS A CEILING OF APPROXIMATELY 
COMMENT	17.5. I HAVE ALLOWED A SMALL MARGIN BELOW THAT IN CASE THE ROUTINE DECLINES 
COMMENT	 IN ACCURACY IN THE NEIGHBORHOOD OF THIS CEILING.
DO IF ((SQRT(#UUC2) > UC2 OR SQRT(#LLC2) > 16.5) AND LC2 > 0) .
+ COMPUTE LC2 = SQRT(#LLC2) .
+ COMPUTE UCDF = #UUCDF .
END IF .
DO IF (SQRT(#UUC2) > UC2 OR UC2 <= TVAL).
+ COMPUTE UC2 = SQRT(#UUC2).
+ COMPUTE LCDF = #LLCDF.
END IF .
EXECUTE .
COMMENT	From T2D.sps at M. J. Smithson’s site.
COMMENT	Edits by K. Wuensch
COMMENT	Two Independent Samples T.
COMPUTE LOWD = LC2/SQRT(N1*N2/(N1+N2)) .
COMPUTE HIGHD = UC2/SQRT(N1*N2/(N1+N2)) .
COMPUTE D = TVAL*SQRT(N1+N2)/SQRT(N1*N2) .
EXECUTE.
"""
    

class SpssWrap(object):
    #This is the SPSSwrap class from spssaux2.py by JKP, IBM
    #     https://www.ibm.com/developerworks/community/files/app#/file/108ce263-b14e-4375-9f70-6d1f733d132b
    #I have included this class directly because spssauz2 is only in SPSS23 and up
    def __init__(self, width=120):
        """width is the approximate maximum in characters for
        a line.  It may be exceeded if an unwrappable item is
        encountered, and in Unicode mode, the width in bytes
        of the lines may be greater than the number of characters.
        However the SPSS 251 limit allows for character expansion
        up to 1000 bytes, which will always be safe.  The default 
        for width is set to facilitate readability, but lines
        containing literals may be chopped up more than the
        width rule requires."""
        
        self.width = width
        self.twrap = textwrap.TextWrapper(width, 
            break_long_words=False, break_on_hyphens=False)
        self.lit = re.compile(r"""('(.*?)'(?!'))|("(.*?)"(?!"))""", re.DOTALL)
        
    def wrap(self, cmdset):
        """Wrap SPSS syntax (not programs) safely and return list of lines
        
        cmdset is the set of commands to wrap as a single string
        or a sequence of cmds.

        cmdset is expected to be an SPSS syntax string
        containing one or more complete commands with line breaks
        separating commands just as used in Submit or a sequence 
        of commands with each command a single item in the list.
        
        Literals will be wrapped as necessary, according to
        SPSS literal rules.  The syntax should not contain literals
        that have already been wrapped with "+"
        
        Unsupported functionality:
            Batch syntax formatting
            GPL code
            Inline comments (/*) other than in the initial position
            Blank line implied command termination
            begin program and begin data blocks"""
        
        lines = []
        # Make each command a separate item in a list of commands
        # Command ends are determined by a trailing period.
        # A command like compute x = 123.\n + 10. will fail, but so
        # does Statistics
        if not spssaux._isseq(cmdset):
            cmdset = cmdset.split(".\n")
            if cmdset[-1] == "":
                del(cmdset[-1])
            for i in range(len(cmdset)):
                if not cmdset[i].endswith("."):
                    cmdset[i] = cmdset[i] + "."
        else:
            # first trim trailing blanks
            for i in range(len(cmdset)):
                cmdset[i] = cmdset[i].rstrip()
            # make sure whole sequence ends with period.
            if not cmdset[-1].endswith("."):
                cmdset[-1] = cmdset[-1] + "."
            # group command fragments in separate lines together
            cmdset2 = []
            fragment = []
            for line in cmdset:
                if line.endswith("."):
                    cmdset2.append("\n".join(fragment) + [line])
                    fragment = []
                else:
                    fragment.append(line)
            cmdset = cmdset2

        for cmd in cmdset:
            loc = 0
            cmdlen = len(cmd)
            if cmdlen > 0:
                # Treat comment commands specially.
                # They continue automatically in most cases, but
                # a /* comment or a wrap point at a period would not continue
                cmdsplit = cmd.split()
                cmdname = cmdsplit[0].lower()
                if cmdname in ["*", "comment", "/*"]:
                    cmdlines = self.twrap.wrap(cmd)
                    for i in range(1, len(cmdlines)):
                        cmdlines[i] = "* " + cmdlines[i]
                    lines.extend(cmdlines)
                    continue
                if cmdname == "begin" and cmdlen > 1 and \
                   cmdsplit[1].lower() in ["data", "program", "gpl"]:
                    raise ValueError("""BEGIN DATA, BEGIN PROGRAM, and BEGIN GPL blocks are not supported in SpssWrap""")
            while loc < cmdlen:
                mo = re.search(self.lit, cmd[loc:])
                if mo is None:
                    lines.extend(self.twrap.wrap(cmd[loc:]))
                    break

                nonlit = cmd[loc:loc + mo.start(0)]
                if len(nonlit) > 0:
                    lines.extend(self.twrap.wrap(nonlit))
                lit = cmd[loc + mo.start(0): loc + mo.end(0)][1:-1]
                litlen = len(lit)
                delim = cmd[loc + mo.start(0)]
                loc += mo.end(0)
                if litlen == 0:
                    lines.append(delim+delim)
                else:
                    for i in range(0, litlen, self.width):
                        fragment = delim + lit[i:i+self.width] + delim
                        if i + self.width < litlen:
                            fragment += "+"
                        lines.append(fragment)
        # clean up any newlines in literals to prevent bad line breaking
        for i in range(len(lines)):
            lines[i] = lines[i].replace("\n", "\\n")
        return lines

def Run(args):
    """Execute the ETTEST extension command"""
    #This is the first function called by the synax command
    #You have to then build a syntax templace and pass it to the processcmd function in SPSS
    #Which will then pass the processed arguments on to the function named do_ettest for the actual
    #Processing
    #print "Entered Run sub of ETTEST"
    #print(args)    
	
    args = args[args.keys()[0]]

    oobj = Syntax([
        Template("GROUPS", subc="",  ktype="existingvarlist", var="groups", islist=True),
        Template("GONE", subc="GROUP", ktype="int", var="gone"),
        Template("GTWO", subc="GROUP", ktype="int", var="gtwo"),       
        Template("DVLIST", subc="WITH", ktype="existingvarlist", var="dvlist", islist=True, vallist=["MISSING","LISTWISE"]),
        Template("MISSING", subc="OPTIONS", ktype="str", var="missing", islist=False),
        Template("CI", subc="OPTIONS", ktype="int", var="conf")
        ])   

    #enable localization
    global _
    try:
        _("---")
    except:
        def _(msg):
            return msg
    # A HELP subcommand overrides all else
    if args.has_key("HELP"):
        #print helptext
        print("Help requested")
    else:
        processcmd(oobj, args, do_ettest, vardict=spssaux.VariableDict())

