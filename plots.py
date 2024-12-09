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


#zad 3 / 24
c3 = TCanvas('cEta', 'Eta', 600, 600)
c3.SetLeftMargin(0.15)  #Space for printing Y label
histo3 = gROOT.FindObject('hEta')
histo3.GetXaxis().SetTitle("Eta")
histo3.GetYaxis().SetTitle("Eta * Charge")
histo3.SetTitle("Muon Pseudorapidity and charge analysis")
fitFunc1=TF1("fitFunc1", "pol1", -2, 2)
histo3.Fit("fitFunc1", "R", " ")
fitFunc1.Draw('same')
histo3.DrawCopy('COL')
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


#zad 17.1
c20 = TCanvas('cVPPGPT1', 'Comparison of Tranverse Momentum from Propagation and Vertex at Station 1 entry', 600, 600)
c20.SetLeftMargin(0.15)  #Space for printing Y label
histo20 = gROOT.FindObject('hVPPGPT1')
histo20.GetXaxis().SetTitle("From Vertex")
histo20.GetYaxis().SetTitle("Propagated")
fitFunc1=TF1("fitFunc1", "pol1", 0, 80)
histo20.Fit("fitFunc1", "R", " ")
fitFunc1.Draw('same')
histo20.DrawCopy('COL')
c20.Print("./plots/"+c20.GetName()+".png")
c20.Update()
input('press enter to exit')

#zad 17.2 
c21= TCanvas('cVPPGPT2', 'Comparison of Tranverse Momentum from Propagation and Vertex at Station 2 entry', 600, 600)
c21.SetLeftMargin(0.15)  #Space for printing Y label
histo21 = gROOT.FindObject('hVPPGPT2')
histo21.GetXaxis().SetTitle("From Vertex")
histo21.GetYaxis().SetTitle("Propagated")
fitFunc1=TF1("fitFunc1", "pol1", 0, 80)
histo21.Fit("fitFunc1", "R", " ")
fitFunc1.Draw('same')
histo21.DrawCopy('COL')
c21.Print("./plots/"+c21.GetName()+".png")
c21.Update()
input('press enter to exit')

#zad 17.3 
c22 = TCanvas('cVSPT1', 'Comparison of Tranverse Momentum from Simulation and Vertex at Station 1 entry', 600, 600)
c22.SetLeftMargin(0.15)  #Space for printing Y label
histo22 = gROOT.FindObject('hVSPT1')
histo22.GetXaxis().SetTitle("From Vertex")
histo22.GetYaxis().SetTitle("Simulated")
fitFunc1=TF1("fitFunc1", "pol1", 0, 80)
histo22.Fit("fitFunc1", "R", " ")
fitFunc1.Draw('same')
histo22.DrawCopy('COL')
c22.Print("./plots/"+c22.GetName()+".png")
c22.Update()
input('press enter to exit')

#zad 17.4 
c23 = TCanvas('cVSPT2', 'Comparison of Tranverse Momentum from Simulation and Vertex at Station 2 entry', 600, 600)
c23.SetLeftMargin(0.15)  #Space for printing Y label
histo23 = gROOT.FindObject('hVSPT2')
histo23.GetXaxis().SetTitle("From Vertex")
histo23.GetYaxis().SetTitle("Simulated")
fitFunc1=TF1("fitFunc1", "pol1", 0, 80)
histo23.Fit("fitFunc1", "R", " ")
fitFunc1.Draw('same')
histo23.DrawCopy('COL')
c23.Print("./plots/"+c23.GetName()+".png")
c23.Update()
input('press enter to exit')

#zad 17.5
c24 = TCanvas('c1Dtest', 'Transverse momentum station 1', 600, 600)
histo24 = gROOT.FindObject('h1Dtest')
histo24.GetXaxis().SetTitle("Simulated/Gen")
#c24.SetLeftMargin(0.15)  #Space for printing Y label
histo24.DrawCopy('COL')
c24.Print("./plots/"+c24.GetName()+".png")
c24.Update()
input('press enter to exit')

#zad 18.0
c26 = TCanvas('c2Dtest', 'Transverse momentum tp.pt 2', 600, 600)
histo26 = gROOT.FindObject('h2Dtest')
histo26.GetXaxis().SetTitle("Simulated/Gen")
#c24.SetLeftMargin(0.15)  #Space for printing Y label
histo26.DrawCopy('COL')
c26.Print("./plots/"+c26.GetName()+".png")
c26.Update()
input('press enter to exit')


#zad 18
c25 = TCanvas('cLandau', 'Landau', 600, 600)
histo25 = gROOT.FindObject('hLandau')
c25.SetLeftMargin(0.15)  #Space for printing Y label
histo25.GetXaxis().SetTitle("Pt_Sim/Pt_Vertex")
histo25.GetYaxis().SetTitle("Energy loss")
histo25.SetTitle("Difference between muon Pt_Sim at station 1 and 2 entries")
fitFuncL=TF1("fitFuncL", "landau", 0, 5)
histo25.Fit("fitFuncL", "R", " ")
#histo25.SetStats(0)
#histo25.GetYaxis().SetRangeUser(0, 1400)
fitFuncL.Draw('same')
histo25.DrawCopy('COL')
c25.Print("./plots/"+c25.GetName()+".png")
c25.Update()
input('press enter to exit')




#zad 19
c27 = TCanvas('cPhiComp', 'Comparison of Phi and PhiB from Hit', 600, 600)
c27.SetLeftMargin(0.15)  #Space for printing Y label
histo27 = gROOT.FindObject('hPhiComp')
histo27.GetXaxis().SetTitle("Phi")
histo27.GetYaxis().SetTitle("PhiB")
#fitFunc1=TF1("fitFunc1", "pol1", -100, 100)
histo27.Fit("fitFunc1", "R", " ")
#fitFunc1.Draw('same')
histo27.DrawCopy('COL')
c27.Print("./plots/"+c27.GetName()+".png")
c27.Update()
input('press enter to exit')



#zad 21
c29 = TCanvas('cPhiBCompSt1', 'PhiB Comb', 600, 600)
c29.SetLeftMargin(0.15)  #Space for printing Y label
histo29 = gROOT.FindObject('hPhiBCompSt1')
histo29.GetYaxis().SetTitle("PhiB Rec")
histo29.GetXaxis().SetTitle("PhiB Sim")
histo29.SetTitle("Comparing PhiB at Station 1")
histo29.DrawCopy('COL')

#fitFunc1=TF1("fitFunc1", "pol1", 0, 300)
#histo29.Fit("fitFunc1", "R", " ")
#fitFunc1.Draw('same')

line = TLine(-1, -1, 1, 1)
line.SetLineColor(2)
#line.SetLineWidth(2)
line.Draw("same")

c29.Print("./plots/"+c29.GetName()+".png")
c29.Update()
input('press enter to exit')

#zad 26
c34 = TCanvas('cPhiBCompSt2', 'PhiB Comb', 600, 600)
c34.SetLeftMargin(0.15)  #Space for printing Y label
histo34 = gROOT.FindObject('hPhiBCompSt2')
histo34.GetYaxis().SetTitle("PhiB Rec")
histo34.GetXaxis().SetTitle("PhiB Sim")
histo34.SetTitle("Comparing PhiB at Station 2")
#fitFunc1=TF1("fitFunc1", "pol1", 0, 300)
#histo29.Fit("fitFunc1", "R", " ")
#fitFunc1.Draw('same')
histo34.DrawCopy('COL')
line = TLine(-1, -1, 1, 1)
line.SetLineColor(2)
#line.SetLineWidth(2)
line.Draw("same")
c34.Print("./plots/"+c34.GetName()+".png")
c34.Update()
input('press enter to exit')


#zad 22
c30 = TCanvas('cHowMany1', 'Phi at St1 Ch2', 600, 600)
c30.SetLeftMargin(0.15)  #Space for printing Y label
histo30 = gROOT.FindObject('hHowMany1')
#histo30.GetXaxis().SetTitle("PhiB Rec")
#histo30.GetYaxis().SetTitle("PhiB Sim")
histo30.SetTitle("Phi at St1 Ch2")
#fitFunc1=TF1("fitFunc1", "pol1", 0, 300)
#histo29.Fit("fitFunc1", "R", " ")
#fitFunc1.Draw('same')
histo30.DrawCopy('COL')
c30.Print("./plots/"+c30.GetName()+".png")
c30.Update()
input('press enter to exit')

#zad 22.5
c31 = TCanvas('cCheck0', 'Comparing vector sizes', 600, 600)
c31.SetLeftMargin(0.15)  #Space for printing Y label
histo31 = gROOT.FindObject('hCheck0')
histo31.GetXaxis().SetTitle("Rec vector")
histo31.GetYaxis().SetTitle("Sim vector")
histo31.SetTitle("Size comparison")
histo31.DrawCopy('COL')
c31.Print("./plots/"+c31.GetName()+".png")
c31.Update()
input('press enter to exit')

#zad 25.1
c32 = TCanvas('cPhiCompareSt1', 'Phi Comp st 1', 600, 600)
c32.SetLeftMargin(0.15)  #Space for printing Y label
histo32 = gROOT.FindObject('hPhiCompSt1')
histo32.GetYaxis().SetTitle("Phi Rec")
histo32.GetXaxis().SetTitle("Phi Sim")
histo32.SetTitle("Comparing Phi at Station 1 entry")
histo32.DrawCopy('COL')
line = TLine(-4, -4, 4, 4)
line.SetLineColor(2)
#line.SetLineWidth(2)
line.Draw("same")
c32.Print("./plots/"+c32.GetName()+".png")
c32.Update()
input('press enter to exit')

#zad 25.2
c33 = TCanvas('cPhiCompareSt2', 'Phi Comp st 2', 600, 600)
c33.SetLeftMargin(0.15)  #Space for printing Y label
histo33 = gROOT.FindObject('hPhiCompSt2')
histo33.GetYaxis().SetTitle("Phi Rec")
histo33.GetXaxis().SetTitle("Phi Sim")
histo33.SetTitle("Comparing Phi at Station 2 entry")
histo33.DrawCopy('COL')
line = TLine(-4, -4, 4, 4)
line.SetLineColor(2)
#line.SetLineWidth(2)
line.Draw("same")
c33.Print("./plots/"+c33.GetName()+".png")
c33.Update()
input('press enter to exit')


#zad 27
c35 = TCanvas('cPhiB_st1', 'PhiB(Pt) St1', 600, 600)
c35.SetLeftMargin(0.15)  #Space for printing Y label
c35.SetLogy(True)
histo35 = gROOT.FindObject('hPhiB_st1')
histo35.GetXaxis().SetTitle("Pt Vertex")
histo35.GetYaxis().SetTitle("PhiB Sim")
histo35.SetTitle("PhiB at station 1 entry")
#fitFunc1=TF1("fitFunc1", "pol1", -100, 100)
#histo28.Fit("fitFunc1", "R", " ")
#fitFunc1.Draw('same')
histo35.DrawCopy('COL')
c35.Print("./plots/"+c35.GetName()+".png")
c35.Update()
input('press enter to exit')


#zad 20
c28 = TCanvas('cPhiB_st2', 'PhiB(Pt) st 2', 600, 600)
c28.SetLeftMargin(0.15)  #Space for printing Y label
c28.SetLogy(True)
histo28 = gROOT.FindObject('hPhiB_st2')
histo28.GetXaxis().SetTitle("Pt Vertex")
histo28.GetYaxis().SetTitle("PhiB Sim")
histo28.SetTitle("PhiB at station 2 entry")
#fitFunc1=TF1("fitFunc1", "pol1", -100, 100)
#histo28.Fit("fitFunc1", "R", " ")
#fitFunc1.Draw('same')
histo28.DrawCopy('COL')
c28.Print("./plots/"+c28.GetName()+".png")
c28.Update()
input('press enter to exit')


#zad 28.1
c36 = TCanvas('cDeltaPhiB1', 'Delta PhiB at station 1 entry', 600, 600)
histo36 = gROOT.FindObject('hDeltaPhiB1')
histo36.SetTitle("Delta PhiB at station 1 entry")
histo36.DrawCopy('COL')
c36.Print("./plots/"+c36.GetName()+".png")
c36.Update()
input('press enter to exit')

#zad 28.2
c37 = TCanvas('cDeltaPhiB2', 'Delta PhiB at station 2 entry', 600, 600)
histo37 = gROOT.FindObject('hDeltaPhiB2')
histo37.SetTitle("Delta PhiB at station 2 entry")
histo37.DrawCopy('COL')
c37.Print("./plots/"+c37.GetName()+".png")
c37.Update()
input('press enter to exit')

#zad 29.1
c38 = TCanvas('cDeltaPhi1', 'Delta Phi at station 1 entry', 600, 600)
histo38 = gROOT.FindObject('hDeltaPhi1')
histo38.SetTitle("Delta Phi at station 1 entry")
histo38.DrawCopy('COL')
c38.Print("./plots/"+c38.GetName()+".png")
c38.Update()
input('press enter to exit')

#zad 29.2
c39 = TCanvas('cDeltaPhi2', 'Delta Phi at station 2 entry', 600, 600)
histo39 = gROOT.FindObject('hDeltaPhi2')
histo39.SetTitle("Delta PhiB at station 2 entry")
histo39.DrawCopy('COL')
c39.Print("./plots/"+c39.GetName()+".png")
c39.Update()
input('press enter to exit')



#zad 30.1 To jest szkic dla konkretnego code 
c40 = TCanvas('cDeltaBCodeSt1', 'DeltaPhiB in the code function', 600, 600)
#c40.SetLeftMargin(0.15)  #Space for printing Y label
histo40 = gROOT.FindObject('hDeltaBCodeSt1')
code_value = 2
bin_code = histo40.GetXaxis().FindBin(code_value)
histo41 = histo40.ProjectionY(f"histo41{code_value}", bin_code, bin_code)
c41 = TCanvas('cSt1Code2', 'Delta PhiB in Station 1 code 2', 600, 600)
histo41.GetXaxis().SetTitle(f'Code value = {code_value}')
histo41.DrawCopy('COL')
c41.Print("./plots/"+c41.GetName()+".png")
c41.Update()
input('press enter to exit')

'''
#zad 30.1
histo40 = gROOT.FindObject('hDeltaBCodeSt1')
for (code_value) in range(2, 7):
    bin_code = histo40.GetXaxis().FindBin(code_value)
    histo_projection = histo40.ProjectionY(f"histo_projection_{code_value}", bin_code, bin_code)

    c40 = TCanvas(f'cSt1Code{code_value}', f'Delta PhiB in Station 1 code {code_value}', 600, 600)
    histo_projection.GetXaxis().SetTitle(f'DeltaPhiB (Code = {code_value})')
    histo_projection.GetYaxis().SetTitle('Entries')
    histo_projection.SetTitle('Delta PhiB in Station 1')
    histo_projection.Draw("COL")

    c40.Print(f"./plots/{c40.GetName()}.png")
    c40.Update()

print("Finished generating projections and saving plots!")

#zad 30.2
histo41 = gROOT.FindObject('hDeltaBCodeSt2')
for (code_value) in range(2, 7):
    bin_code = histo41.GetXaxis().FindBin(code_value)
    histo_projection = histo41.ProjectionY(f"histo_projection_{code_value}", bin_code, bin_code)

    c41 = TCanvas(f'cSt2Code{code_value}', f'Delta PhiB in Station 2 code {code_value}', 600, 600)
    histo_projection.GetXaxis().SetTitle(f'DeltaPhiB (Code = {code_value})')
    histo_projection.GetYaxis().SetTitle('Entries')
    histo_projection.SetTitle('Delta PhiB in Station 2')
    histo_projection.Draw("COL")

    c41.Print(f"./plots/{c41.GetName()}.png")
    c41.Update()

print("Finished generating projections and saving plots!")
