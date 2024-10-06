#!/cvmfs/cms.cern.ch/el8_amd64_gcc11/cms/cmssw/CMSSW_13_2_8/external/el8_amd64_gcc11/bin/python3

#Original code of plots.py
import sys
import math
from ROOT import *



#print ("Hello ROOT")
fileName = "licDigiHistos.root"

print ('Read data from: ', fileName)
gROOT.Reset()
f = TFile(fileName) #;
f.ls()  #;

'''
#przyklad
c0 = TCanvas('cLicExample','cLicExample',600,600)
histo0 = gROOT.FindObject('hLicExample')
histo0.DrawCopy('box text')
c0.Print("./plots/"+c0.GetName()+".png")
c0.Update()
input('press enter to exit')


#zad 1
c1 = TCanvas('cPt','Histogram of Transverse Momentum',600,600)
histo1 = gROOT.FindObject('hPt')
histo1.DrawCopy('hist text')
c1.Print("./plots/"+c1.GetName()+".png")
c1.Update()
input('press enter to exit')


#zad 2
c2 = TCanvas('cVz','Histogram of Z0-Vertex',600,600)
histo2 = gROOT.FindObject('hVz')
fitFunc=TF1("fitFunc","gaus", -15, 15) #Fitting Function #zad 6 fitting 
histo2.Fit("fitFunc", "R", " ")
histo2.DrawCopy('hist text') 
fitFunc.Draw("same")
c2.Print("./plots/"+c2.GetName()+".png")
c2.Update()
input('press enter to exit')


#zad 3
c3 = TCanvas('cEta','Histogram of Eta',600,600)
histo3 = gROOT.FindObject('hEta')
histo3.DrawCopy('hist text')
c3.Print("./plots/"+c3.GetName()+".png")
c3.Update()
input('press enter to exit')



#zad 4
c4 = TCanvas('cVxy', 'Histogram of Vertex X and Y', 600, 600)
histo4 = gROOT.FindObject('hVxy')
histo4.DrawCopy('box text')
c4.Print("./plots/"+c4.GetName()+".png")
c4.Update()
input('press enter to exit')

#zad 5
c5 = TCanvas('cT', 'Determining time unit', 600, 600)
histo5 = gROOT.FindObject('hT')
fitFunc1=TF1("fitFunc1", "pol1", 0, 10e-8)
histo5.Fit("fitFunc1", "R", " ")
#The following range manipulation wasn't proven successful
# histo5.GetXaxis().SetRangeUser(0, 6*10e-8) # Set range of x-axis 
# histo5.GetYaxis().SetRangeUser(0, 6*10e-8) # Set range of y-axis
fitFunc1.Draw("same")
histo5.DrawCopy('box text')
c5.Print("./plots/"+c5.GetName()+".png")
c5.Update()
input('press enter to exit')



#zad 7.1
c6 = TCanvas('cPSimHit', 'Histogram of numberOfHits', 600, 600)
histo6 = gROOT.FindObject('hPSimHit')
histo6.GetXaxis().SetRangeUser(0, 80) 
histo6.DrawCopy('hist text')
c6.Print("./plots/"+c6.GetName()+".png")
c6.Update()
input('press enter to exit')

#zad 7.2
c7 = TCanvas('cPSimTrackerHit', 'Histogram of numberOfTrackerHits', 600, 600)
histo7 = gROOT.FindObject('hPSimTrackerHit')
histo7.GetXaxis().SetRangeUser(0, 30) 
histo7.DrawCopy('hist text')
c7.Print("./plots/"+c7.GetName()+".png")
c7.Update()
input('press enter to exit')

#zad 7.3
c8 = TCanvas('cPSimHitVector', 'Histogram of PSimHit Size from vector', 600, 600)
histo8 = gROOT.FindObject('hPSimHitVector')
#histo8.GetXaxis().SetRangeUser(0, 30) 
histo8.DrawCopy('hist text')
c8.Print("./plots/"+c8.GetName()+".png")
c8.Update()
input('press enter to exit')

#zad 8
c9 = TCanvas('cPSimHitXYZ', 'Histogram of PSimHit Position', 600, 600)
histo9 = gROOT.FindObject('hPSimHitXYZ')
#histo9.DrawCopy('BOX')
histo9.DrawCopy('LEGO')
c9.Print("./plots/"+c9.GetName()+".png")
c9.Update()
input('press enter to exit')


#zad 11
c10 = TCanvas('cPPGvsS', 'Histogram of differences between propagation and simulation', 600, 600)
histo10 = gROOT.FindObject('hPPGvS')
histo10.DrawCopy('hist text')
c10.Print("./plots/"+c10.GetName()+".png")
c10.Update()
input('press enter to exit')

#zad 13
c11 = TCanvas('cPSimHitRZ', 'PSimHit Map', 600, 600)
histo11 = gROOT.FindObject('hPSimHitRZ')
histo11.GetXaxis().SetTitle("R")
histo11.GetYaxis().SetTitle("Z")
histo11.DrawCopy('box text')
c11.Print("./plots/"+c11.GetName()+".png")
c11.Update()
input('press enter to exit')

#zad 16
c12 = TCanvas('cPSimHitXY', 'PSimHit Map', 600, 600)
histo12 = gROOT.FindObject('hPSimHitRZ')
histo12.GetXaxis().SetTitle("X")
histo12.GetYaxis().SetTitle("Y")
histo12.GetXaxis().SetRangeUser(-800, 800) 
histo12.GetYaxis().SetRangeUser(-800, 800) 
histo12.DrawCopy('box text')
c12.Print("./plots/"+c12.GetName()+".png")
c12.Update()
input('press enter to exit')



#zad 17
c13 = TCanvas('cPHPT', 'Comparison of Tranverse Momentum from Propagation and Simulation at Chambers', 600, 600)
c13.SetLeftMargin(0.15)  #Space for printing Y label
histo13 = gROOT.FindObject('hPHPT')
histo13.GetXaxis().SetTitle("Propagated")
histo13.GetYaxis().SetTitle("Simulated")
fitFunc2=TF1("fitFunc1", "pol1", -20, 20)
histo13.Fit("fitFunc1", "R", " ")
fitFunc2.Draw('same')
histo13.DrawCopy('COL')
c13.Print("./plots/"+c13.GetName()+".png")
c13.Update()
input('press enter to exit')


#zad 18.1
c14 = TCanvas('cPvSX1', 'The difference in X plane position obtained from propagation and simulation', 600, 600)
histo14 = gROOT.FindObject('hPvSX1')
histo14.DrawCopy('COL')
c14.Print("./plots/"+c14.GetName()+".png")
c14.Update()
input('press enter to exit')

#zad 18.2
c15 = TCanvas('cPvSX2', 'The difference in X plane position obtained from propagation and simulation', 600, 600)
histo15 = gROOT.FindObject('hPvSX2')
histo15.DrawCopy('COL')
c15.Print("./plots/"+c15.GetName()+".png")
c15.Update()
input('press enter to exit')

#zad 18.3
c16 = TCanvas('cPvSY1', 'The difference in Y plane position obtained from propagation and simulation', 600, 600)
histo16 = gROOT.FindObject('hPvSY1')
histo16.DrawCopy('COL')
c16.Print("./plots/"+c16.GetName()+".png")
c16.Update()
input('press enter to exit')

#zad 18.4
c17 = TCanvas('cPvSY2', 'The difference in Y plane position obtained from propagation and simulation', 600, 600)
histo17 = gROOT.FindObject('hPvSY2')
histo17.DrawCopy('COL')
c17.Print("./plots/"+c17.GetName()+".png")
c17.Update()
input('press enter to exit')
'''

#zad 19
c18 = TCanvas('cPtX', 'Full width at half maximum of X histogram in relation to Pt value', 600, 600)
histo18 = gROOT.FindObject('hPtX')
histo18.DrawCopy('E')
c18.Print("./plots/"+c18.GetName()+".png")
c18.Update()
input('press enter to exit')

#zad 20
c19 = TCanvas('cSPt', 'Spread in relation to Pt', 600, 600)
histo19 = gROOT.FindObject('hSPt')
histo19.DrawCopy('COL')
c19.Print("./plots/"+c19.GetName()+".png")
c19.Update()
input('press enter to exit')

