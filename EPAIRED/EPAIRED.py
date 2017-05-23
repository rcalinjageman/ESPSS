# -*- coding: utf-8 -*-

"""EPAIRED extension command"""

__author__ =  'RC-J'
__version__=  '1.0.0'

import spss, spssaux, spssdata
from extension import Template, Syntax, processcmd
import random
from xml.etree import ElementTree
import re
import math


# This it the python module for an extension that completes a paired t-test with an estimation approach (EPAIRED)
# The improvements are: reporting cohen's d and the CI for cohen's d.
# Cohen's d is reported correct for bias (dunbiased = Hedges g)

# This module was written by Bob Calin-Jageman

# This module contains code adapted from:
#       ESCI by Geoff Cumming
#       Python documentation
#       And lots of inspiration from the extension STATS_CORRELATION.py

# To do:
#	Check on why Pearson's r caclulation is a bit off; consider slurping it directly from the SPSS results
#       Check what happens when N > 200
#       Add a figure
#       Integrate into one package with ETTEST and CORRELATION

def do_epaired(plist, wlist, missing, conf):
    #This is the actual function that completes the EPAIRED procedure
      #It accepts the list of arguments parsed out in the run subroutine SPSS initially calls
      #The run subroutine calls this function via SPSS' processcmd routine

    #Set to True to get debugging output; set to False for distribution
    debugit = False
    
    if debugit:
        print("Entered do_paired")
        print(plist)
        print(wlist)
        print missing
        print conf

    #Get the name for the current dataset
    activeds = spss.ActiveDataset()

    #Pre-req checks - Make sure all requirements are met to run the command
    #To do: write this section


    #Put the parameters for the syntax command back into string format
    paired_string = " ".join(plist)
    with_string = " ".join(wlist)

    #Setup the t-test command--will need to adapt this for the more complex t-test commands possible
    mycommand="""T-TEST PAIRS=%(paired_string)s WITH %(with_string)s (PAIRED)
/MISSING=ANALYSIS
/CRITERIA=CI(0.%(conf)s).""" % locals()  
    if debugit:
        print(mycommand)

    #Create the descriptive stats output
    tdestable,err = spssaux.CreateDatasetOutput(mycommand, omsid='T-Test', subtype='Paired Samples Statistics')
    #Probably should do some checking on the spss failcode as well    
    if not err == 0:
        raise ValueError(_("Descriptive Stats Command Failed"))
    spss.Submit("DATASET ACTIVATE "+tdestable+".")

    #Now compute t crit for each individual mean so we can gather CIs
    ci_level = float(conf)/100
    ci_level = ci_level + ((1-ci_level)/2)
    
    gettcrit = """COMPUTE mtcrit = IDF.T("""+str(ci_level)+""", N-1).
EXECUTE."""
    if debugit:
        print(gettcrit)
    spss.Submit(gettcrit)

    desdata = spssdata.Spssdata()
    resultlist = []
    start = True
    count = 0

    for row in desdata:
        if start:
            resultlist.append(pairedt())
            resultlist[count].conf = conf
            resultlist[count].m1 = row.Mean
            resultlist[count].n = row.N
            resultlist[count].s1 = row[7]
            resultlist[count].sem1 = row[8]
            resultlist[count].m1_tcrit = row.mtcrit
            resultlist[count].var1 = row.Var2
            start = False
        else:
            resultlist[count].m2 = row.Mean
            resultlist[count].s2 = row[7]
            resultlist[count].sem2 = row[8]
            resultlist[count].m2_tcrit = row.mtcrit
            resultlist[count].var2 = row.Var2
            start = True
            count = count + 1

    desdata.CClose()
    spss.Submit("DATASET ACTIVATE "+activeds+".")
    spss.Submit("DATASET CLOSE "+tdestable+".")

    #Now re-do the command to collect the t-table output
    ttable,err = spssaux.CreateDatasetOutput(mycommand, omsid='T-Test', subtype='Paired Samples Test')
    #Probably should do some checking on the spss failcode as well    
    if not err == 0:
        raise ValueError(_("T-TEST Command Failed"))
    spss.Submit("DATASET ACTIVATE "+ttable+".")

    desdata = spssdata.Spssdata()
    count = 0

    for row in desdata:
        #Remains stupid that some columns can't be addressed by name due to illeagal characters.  
        resultlist[count].mdiff = row.Mean
        resultlist[count].sdiff = row[6]
        resultlist[count].mdiff_sem = row[7]
        resultlist[count].mdiff_low = row.Lower
        resultlist[count].mdiff_high = row.Upper
        resultlist[count].t = row.t
        resultlist[count].df = row.df
        resultlist[count].p = row[12]
        count = count + 1

    desdata.CClose()
    spss.Submit("DATASET ACTIVATE "+activeds+".")
    spss.Submit("DATASET CLOSE "+ttable+".")

    #Don't forget to patch up the results from SPSS, filling in missing info (CI on d, MoE and CI for individual means, MoE for MDiff
    for myresult in resultlist:
            myresult.calc_missing()
            if debugit:
                myresult.printout()
    display(resultlist)
    
def display(resultlist):
    #Do the actual output of the epaired
    spss.StartProcedure(_("Estimation: Compare two related means"), "EPAIRED")

    for myresult in resultlist:
        myrowlabels = [myresult.var1, myresult.var2, "Raw Score Difference", "Standardized difference (Cohen's d)"]
        mycolumnlabels = ["M", str(myresult.conf)+"% MoE", str(myresult.conf)+"% CI [low, high]", "s"]

        m1ci = "[{0:.3f}, {1:.3f}]".format(myresult.m1_low, myresult.m1_high)
        m2ci = "[{0:.3f}, {1:.3f}]".format(myresult.m2_low, myresult.m2_high)
        dci = "[{0:.3f}, {1:.3f}]".format(myresult.d_low, myresult.d_high)
        mci = "[{0:.3f}, {1:.3f}]".format(myresult.mdiff_low, myresult.mdiff_high)
        if myresult.conf != 95:
            dci = "95% CI "+dci
        
        mycells = [[myresult.m1, myresult.m1_moe, m1ci, myresult.s1],
                     [myresult.m2, myresult.m2_moe, m2ci, myresult.s2],
                     [myresult.mdiff, myresult.mdiff_moe, mci, "-"],
                     [myresult.d, myresult.d_moeav, dci, "-"]
                     ]
        mycaption = "N = {0:.0f}; correlation between measures = {1:.3f}".format(myresult.n, myresult.r)
        
        mytable = spss.BasePivotTable("Estimate difference from " + myresult.var1 + " to " + myresult.var2, "EPAIRED", caption=mycaption)
        mytable.SimplePivotTable(rowdim = "Row",
                                      rowlabels= myrowlabels,
                                      coldim = "Column",
                                      collabels = mycolumnlabels,
                                      cells = mycells
                                      )
        eq_tresult =  "NHST Approach: t({0:.0f}) = {1:.2f}, p = {2:.4f}.".format(myresult.df, myresult.t, myresult.p)
        eq_txtblock = spss.TextBlock("NHST Approach comparing "+ myresult.var1 + " to " + myresult.var2, eq_tresult)

    othercomments = spss.TextBlock("Important Notes", """

The EPAIRED routine was written by Bob Calin-Jageman, v1.0, 5/23/2017
* Cohen's d is calculated with correction for bias (d-unbiased).  This is also often referred to as Hedge's g.
* For Cohen's d, the CI is not symmetrical.  The MoE reported is 1/2 the distance between the upper and lower bounds (MoE-average)
* For Cohen's d, the CI is **ALWAYS** with a 95% confidence level, even if a different confidence interval was requested.  This is because the technique for calculate the CI of d is only vetted for a 95% confidence level.
* The CI for Cohen's d is calculated using the method of Algina & Kesselman (2003)
* The calculation for Pearson's r seems to be just a bit off
* The method for calculating CI for Cohen's d may break down when N > 200... need to do some further testing on this.

Some tips:
* Visualize your data--don't draw conclusions until you've made a good visualization showing
all the data that emphasizes the group comparison.  Try ESCI for this.
* Don't draw causal conclusions unless counter-balancing was used.
* The confidence intervals only make sense to the degree to which the sample
is representative of the population.
* Clearly predicted (confirmatory) results are to be trusted more than those found by
exploring the data, but in either case seek replication.
* Be sure to present the whole story when reporting results; statistics can only be judged in
the context of all the data collected and analyses run.
""")
        
    spss.EndProcedure()

class pairedt:
        conf = 0
        n = 0
        m1 = 0
        m2 = 0
        m1_tcrit = 0
        m2_trcit = 0
        m1_low = 0
        m2_low = 0
        m1_high = 0
        m2_high = 0
        m1_moe = 0
        m2_moe = 0
        s1 = 0
        s2 = 0
        sem1 = 0
        sem2 = 0
        d = 0
        r = 0
        d_low = 0
        d_high = 0
        d_moeav = 0
        mdiff = 0
        mdiff_sem = 0
        mdiff_moe = 0
        mdiff_low = 0
        mdiff_high = 0
        sdiff = 0
        var1 = "Var1"
        var2 = "Var2"
        t = 0
        p = 0
        df = 0
        flagged = False
        
        def printout(self):
            print "Compare paired means report"
            print "Var1 = ", self.var1, " Var2 = ", self.var2
            for property, value in vars(self).iteritems():
                print property, ": ", value
            print "d = ", self.d, " 95%CI [", self.d_low, ", ", self.d_high, "]"
            print "End Report"
            print

        def calc_missing(self):
            self.m1_moe = self.m1_tcrit * self.sem1
            self.m2_moe = self.m2_tcrit * self.sem2
            self.m1_low = self.m1 - self.m1_moe
            self.m2_low = self.m2 - self.m2_moe
            self.m1_high = self.m1 + self.m1_moe
            self.m2_high = self.m2 + self.m2_moe
            self.mdiff_moe = (self.mdiff_high - self.mdiff_low) / 2
            myresult = cipaired(self.m1, self.m2, self.s1, self.s2, self.n, self.sdiff, 0.95)
            self.d = myresult[1]
            self.d_low = myresult[2]
            self.d_high= myresult[3]
            self.d_moeav = (self.d_high - self.d_low) / 2
            self.flagged = myresult[4]
            dbiased = myresult[0]
            self.r = 1 - ((dbiased * dbiased * self.n) / (2 * self.t * self.t))
            self.var1 = self.var1.strip()
            self.var2 = self.var2.strip()
            
    

def nct(tval, delta, df):
  #This is the non-central t-distruction function
  #This is a direct translation from visual-basic to Python of the nct implementation by Geoff Cumming for ESCI

  #Noncentral t function, after Excel worksheet calculation by Geoff Robinson
  #The three parameters are:
  #   tval -- our t value
  #   delta -- noncentrality parameter
  #   df -- degrees of freedom

  #Dim dfmo As Integer     'for df-1
  #Dim tosdf As Double     'for tval over square root of df
  #Dim sep As Double       'for separation of points
  #Dim sepa As Double      'for temp calc of points
  #Dim consta As Double    'for constant
  #Dim ctr As Integer      'loop counter

  dfmo = df - 1.0
  tosdf = tval / math.sqrt(df)
  sep = (math.sqrt(df) + 7) / 100
  consta = math.exp( (2 - df) * 0.5 * math.log1p(1) - math.lgamma(df/2) ) * sep /3

  #now do first term in cross product summation, with df=0 special
  if dfmo > 0:
    nctvalue = stdnormdist(0 * tosdf - delta) * 0 ** dfmo * math.exp(-0.5* 0 * 0)
  else:  
    nctvalue = stdnormdist(0 * tosdf - delta) * math.exp(-0.5 * 0 * 0)

  #now add in second term, with multiplier 4
  nctvalue = nctvalue + 4 * stdnormdist(sep*tosdf - delta) * sep ** dfmo * math.exp(-0.5*sep*sep)    

  #now loop 49 times to add 98 terms
  for ctr in range(1,49):
      sepa = 2 * ctr * sep
      nctvalue = nctvalue + 2 * stdnormdist(sepa * tosdf - delta) * sepa ** dfmo * math.exp(-0.5 * sepa * sepa)
      sepa = sepa + sep
      nctvalue = nctvalue + 4 * stdnormdist(sepa * tosdf - delta) * sepa ** dfmo * math.exp(-0.5 * sepa * sepa)

  #add in last term
  sepa = sepa + sep
  nctvalue = nctvalue + stdnormdist(sepa * tosdf - delta) * sepa ** dfmo * math.exp(-0.5 * sepa * sepa)

  #multiply by the constant and we've finished
  nctvalue = nctvalue * consta
  
  return nctvalue

def stdnormdist(x):
  #Cumulative standard normal distribution function
  #Lifted from https://docs.python.org/3.2/library/math.html
  return (1.0 + math.erf(x / math.sqrt(2.0))) / 2.0

def cipaired(m1, m2, s1, s2, n, sdiff, cl):
    tolerance = 0.000001
    maxcycles = 10000
    
    utarget = (1 - cl) / 2
    ltarget = 1- utarget
    df = n - 1.0
    mdiff = m1 - m2
    estsigma = math.sqrt(((s1**2) + (s2**2) ) / 2)
    cohend = mdiff / estsigma
    sem = sdiff / math.sqrt(n)
    pairedt = mdiff / sem
    flagmaxed = False

    #print "mdiff = ", mdiff
    #print "estsigma = ", estsigma
    #print "sem = ", sem
    #print "df = ", df
    #print "paired t = ", pairedt

    unbiasfactor = math.exp( math.lgamma(df / 2 )) / ( math.sqrt(df/2) * math.exp(math.lgamma(df/2 - 0.5)) )
    cohend_unbiased = cohend * unbiasfactor

    approx_low_raw = mdiff - (2.2 * sem)
    approx_high_raw = mdiff + (2.2 * sem)
    approx_low_d = approx_low_raw / estsigma
    approx_high_d = approx_high_raw / estsigma

    low_nct = approx_low_d * estsigma * math.sqrt(n) / sdiff
    high_nct = approx_high_d * estsigma * math.sqrt(n) / sdiff

    #print "starting low_nct = ", low_nct
    #print "starting high_nct = ", high_nct

    cycles = 0
    low_nct_p = nct(pairedt, low_nct, df)
    while (math.fabs(ltarget - low_nct_p) > tolerance) & (cycles < maxcycles):
      #print "Current difference", (low_nct_p - ltarget)
      change = low_nct_p - ltarget
      #print "low_nct = ", low_nct, "; low_nct_p = ", low_nct_p, "; change = ", change
      low_nct = low_nct + (change * math.fabs(low_nct))
      low_nct_p = nct(pairedt, low_nct, df)
      cycles = cycles +1

    if (cycles >= maxcycles):
      flagmaxed = True

    #print "cycles = ", cycles

    high_nct_p = nct(pairedt, high_nct, df)
    cycles = 0
    while (math.fabs(utarget - high_nct_p) > tolerance) & (cycles < maxcycles):
      #print "Current difference", (high_nct_p - utarget)
      change = (high_nct_p - utarget)
      #print "high_nct = ", high_nct, "; high_nct_p = ", high_nct_p, "; change = ", change
      high_nct = high_nct + (change * math.fabs(high_nct))
      high_nct_p = nct(pairedt, high_nct, df)
      cycles = cycles + 1

    if (cycles >= maxcycles):
      flagmaxed = True

    #print "cycles = ", cycles 


    final_low_d = low_nct * sdiff / (estsigma * math.sqrt(n) )
    final_high_d = high_nct * sdiff / (estsigma * math.sqrt(n) )

    #print "approx_low_raw = ", approx_low_raw
    #print "approx_high_raw = ", approx_high_raw
    #print "approx_low_d = ", approx_low_d
    #print "approx_high_d = ", approx_high_d
    #print "low_nct = ", low_nct
    #print "high_nct = ", high_nct
    #print "d = ", cohend
    #print "unbiasfactor = ", unbiasfactor
    #print "cohend_unbiased = ", cohend_unbiased
    #print "final_low_d = ", final_low_d
    #print "final_high_d = ", final_high_d

    return [cohend, cohend_unbiased, final_low_d, final_high_d, flagmaxed]            

def Run(args):
    """Execute the EPAIRED extension command"""
    #This is the first function called by the synax command
    #You have to then build a syntax templace and pass it to the processcmd function in SPSS
    #Which will then pass the processed arguments on to the function named do_ettest for the actual
    #Processing
    #print "Entered Run sub of EPAIRED"
    #print(args)    
	
    args = args[args.keys()[0]]

    oobj = Syntax([
	Template("PAIRS", subc="", ktype="existingvarlist", var="plist", islist=True, ),
        Template("WITH", subc="", ktype="existingvarlist", var="wlist", islist=True, ),
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
        processcmd(oobj, args, do_epaired, vardict=spssaux.VariableDict())

