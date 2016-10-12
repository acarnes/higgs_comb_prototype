from ROOT import *
 
f0 = TFile.Open('simulated_mass_distribution.root')
f1 = TFile('test.root','RECREATE')
 
h_sigVBF = f0.Get('histos/dimuMass_VBF').Clone('sigVBF')
h_sigGGF = f0.Get('histos/dimuMass_GGF').Clone('sigGGF')
 
h_bkgTTJ = f0.Get('histos/dimuMass_TTJets').Clone('bkgTTJ')
h_bkgDYJ = f0.Get('histos/dimuMass_DYJetsToLL').Clone('bkgDYJ')
 
h_data   = f0.Get('Mass').Clone('data_obs')
 
f1.cd()
 
h_sigVBF.Write()
h_sigGGF.Write()
h_bkgTTJ.Write()
h_bkgDYJ.Write()
h_data.Write()
 
f1.Write()
f1.Close()
