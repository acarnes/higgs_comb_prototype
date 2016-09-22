#!/usr/bin/env python
#------------------------------------------------------------------
# Prototype limit setting - Fit simple model to CMS Run2 H to Mu Mu
#------------------------------------------------------------------
import os,sys
from time import sleep
#from histutil import *
from ROOT import *
#------------------------------------------------------------------
def main():
    # set up standard graphics style (see python/histutil.py)
    #setStyle()
    gStyle.SetTitleYOffset(1.75)
    
    #----------------------------------------
    # Read in histogram
    #----------------------------------------
    f = TFile('simulated_mass_distribution.root')
    h = f.Get('Mass')
   
    nbins   = h.GetNbinsX()
    massmin = h.GetBinLowEdge(1)
    massmax = massmin + nbins*h.GetBinWidth(1)

    #----------------------------------------
    # create a workspace so that we can use
    # its factory method
    #----------------------------------------
    # suppress all messages except those that matter
    RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
    print "="*78    
    wspace = RooWorkspace('hmumu')

    #----------------------------------------
    # create di-muon mass variable
    # syntax:
    # <name>[initial-val, min-val, max-val]
    #----------------------------------------
    wspace.factory('x[120.0, %f, %f]' % (massmin, massmax))
    wspace.var('x').SetTitle('m_{#mu#mu}')
    wspace.var('x').setUnit('GeV')

    # define the set obs = (x)
    wspace.defineSet('obs', 'x')

    # make the set obs known to Python
    obs  = wspace.set('obs')

    # create binned dataset from histograms.
    # each histogram in hmap is associated with a diphoton
    # category.
    data = RooDataHist('data_obs', 'data_obs', RooArgList(obs), h)
    
    # add data to workspace.
    # the RooCmdArg() is a workaround a PyROOT "feature"
    getattr(wspace,'import')(data, RooCmdArg())

    #----------------------------------------
    # create background model
    #----------------------------------------
    ndata = int(data.sumEntries());
    wspace.factory('bmodel_norm[%d, 0.0, %d]' % (ndata, 3.0*ndata))
    wspace.factory('a1[ 5.0, -50, 50]')
    wspace.factory('a2[-1.0, -50, 50]')

    # exp(-(a1*(x/100) + a2*(x/100)^2))
    wspace.factory('expr::f("-(a1*(x/100)+a2*(x/100)^2)",a1,a2,x)')
    # exp(c*x), c = 1 and x = f
    wspace.factory('Exponential::bmodel(f, 1)')
    bmodel  = wspace.pdf('bmodel')

    #----------------------------------------
    # create signal model
    #----------------------------------------
    wspace.factory('smodel_norm[100, 0.0, 1000.0]')
    wspace.factory('mass[125, %f, %f]' % (massmin, massmax))
    wspace.factory('w[1.0, 0.1, 10]')
    wspace.factory('Gaussian::smodel(x, mass, w)')
    wspace.var('mass').setConstant()

    #----------------------------------------
    # save data and signal & bg models for use
    # with higgs combine
    #----------------------------------------
    wspace.SaveAs("p.root")

    #----------------------------------------
    # create background + signal model
    # use this to fit the pseudodata and make a picture
    #----------------------------------------
    # According to Wouter, the following automatically
    # produces an extended pdf 
    # p(x|...) = exp(-(background+signal)) * PROD p(x_i|...)
    wspace.factory('SUM::model(bmodel_norm*bmodel, smodel_norm*smodel)')
    model = wspace.pdf('model')

    #----------------------------------------
    # Fit
    #----------------------------------------
    print "fit model to data"
    print "="*80

    b = data.sumEntries()
    wspace.var('bmodel_norm').setVal(b)
    
    swatch = TStopwatch()
    swatch.Start()
    model.fitTo(data)
    print "real time: %10.3f s" % swatch.RealTime()
    
    vbkg = wspace.var('bmodel_norm')
    vsig = wspace.var('smodel_norm')
    vmass= wspace.var('mass')
    vwidth=wspace.var('w')
    
    # Get results and compute a simple measure of signal significance
    bkg   = vbkg.getVal()
    ebkg  = vbkg.getError()
    sig   = vsig.getVal()
    esig  = vsig.getError()
    mass  = vmass.getVal()
    emass = vmass.getError()
    width = vwidth.getVal()
    ewidth= vwidth.getError()
    zvalue= sig / esig

    print "="*80
    print "background: %10.1f +\-%-5.1f GeV" % (bkg, ebkg)
    print "signal:     %10.1f +\-%-5.1f" % (sig, esig)
    print "mass:       %10.1f +\-%-4.1f GeV" % (mass, emass)
    print "width:      %10.1f +\-%-4.1f GeV" % (width, ewidth)
    print "sig/esig:   %10.1f" % zvalue
    print

    wspace.Print()

    x = wspace.var('x')
    xframe = x.frame()
    xframe.GetXaxis().SetNdivisions(505)
    data.plotOn(xframe)
    model.plotOn(xframe)
    bmodel.plotOn(xframe, RooFit.LineColor(kRed+1))

    # hack to prevent display of background parameters
    wspace.var('a1').setConstant()
    wspace.var('a2').setConstant()
    model.paramOn(xframe)

    c1 = TCanvas('fig_hmumu_fit', 'fit', 10, 10, 500, 500)
    xframe.Draw()
    c1.SaveAs('.png')
    
    sleep(5)

#    #-----------------------------------------------------
#    # Define a few useful sets.
#    #-----------------------------------------------------
#    sets = [('obs',  'N1,N2,N3'),    # observations
#            # parameter of interest
#            ('poi',  'mu'),
#            # nuisance parameters
#            ('nuis', 's1,s2,s3,b_zz1,b_zz2,b_zz3,b_zx1,b_zx2,b_zx3')]
#    for t in sets:
#        name, parlist = t
#        wspace.defineSet(name, parlist)

#    #-----------------------------------------------------
#    # Create model configuration. This is needed for the
#    # statistical analyses
#    #-----------------------------------------------------
#    cfg = RooStats.ModelConfig('cfg')
#    cfg.SetWorkspace(wspace)
#    cfg.SetPdf(wspace.pdf('model'))
#    cfg.SetParametersOfInterest(wspace.set('poi'))
#    cfg.SetNuisanceParameters(wspace.set('nuis'))
#
#    # import model configuration into workspace
#    getattr(wspace, 'import')(cfg)

#    #-----------------------------------------------------    
#    # Fit model to data
#    #-----------------------------------------------------
#    results = wspace.pdf('model').fitTo(data, RooFit.Save())
#    results.Print()
#    # get maximum likelihood estimate (MLE) of mu
#    bestFit = wspace.var('mu').getVal()
#
#    #-----------------------------------------------------    
#    # Compute interval based on profile likelihood
#    #-----------------------------------------------------
#    # suppress some (apparently) innocuous warnings
#    msgservice = RooMsgService.instance()
#    msgservice.setGlobalKillBelow(RooFit.FATAL)
#
#    print '== Profile likelihood calculations =='
#    plc = RooStats.ProfileLikelihoodCalculator(data, cfg)
#    CL  = 0.683
#    plc.SetConfidenceLevel(CL)
#    plcInterval = plc.GetInterval()
#    # plot interval
#    plcplot = RooStats.LikelihoodIntervalPlot(plcInterval)
#    plccanvas = TCanvas('fig_PL', 'plc', 10, 10, 400, 400)
#    plcplot.Draw()
#    plccanvas.Update()
#    plccanvas.SaveAs('.png')

#    # get central interval
#    lowerLimit1  = plcInterval.LowerLimit(wspace.var('mu'))
#    upperLimit1  = plcInterval.UpperLimit(wspace.var('mu'))
#    print '\t%4.1f%s CL interval = [%5.3f, %5.3f]' % \
#      (100*CL, '%', lowerLimit1, upperLimit1)
#
#    # compute a 95% upper limit on mu by
#    # computing a 90% central interval and
#    # ignoring the lower limit
#    CL = 0.90
#    plc.SetConfidenceLevel(CL)
#    plcInterval = plc.GetInterval()
#    upperLimit  = plcInterval.UpperLimit(wspace.var('mu'))
#    CL = 0.95
#    print '\t%4.1f%s upper limit = %5.3f' % \
#      (100*CL, '%', upperLimit)
#
#    # compute upper and lower interval bounds relative to
#    # best fit value of mu
#    upper   = upperLimit1-bestFit
#    lower   = bestFit-lowerLimit1
#    print '\tmu (MLE)          = %5.3f -%-5.3f/+%-5.3f' % \
#      (bestFit, lower, upper)

#    # compute a standard measure of significance:
#    #   q(mu) = -2*log [Lp(mu) / Lp(mu_hat)]
#    # where Lp(mu) is the profile likelihood and mu_hat
#    # is the MLE of mu (the best fit value).
#    # basically, we are testing the hypothesis that
#    # mu = 0. We define the significance measure so
#    # that large values cast doubt on the mu=0
#    # hypothesis.
#    # 1. get the negative log-profilelikelihood ratio
#    logplratio = plcInterval.GetLikelihoodRatio()
#    # 2. set value of mu to zero (no-signal)
#    wspace.var('mu').setVal(0)
#    # 3. compute q(0) 
#    q0 = 2*logplratio.getVal()
#    # 4. compute a measure that is approximately Z
#    # standard deviations from mu = 0
#    Z  = sqrt(q0)
#    print "\tZ = sqrt(q(0))    = %5.3f" % Z
#    print "\twhere q(mu) = -2*log[Lp(mu)/Lp(mu_hat)]\n"


#------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print
    print "bye cruel world!"
    print



