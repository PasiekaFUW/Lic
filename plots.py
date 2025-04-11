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


#zad 27
c35 = TCanvas('cPhiB_st1', 'PhiB(Pt) St1', 600, 600)
c35.SetLeftMargin(0.15)  #Space for printing Y label
c35.SetLogy(True)
histo35 = gROOT.FindObject('hPhiB_st1')
histo35.GetXaxis().SetTitle("Transverse Momentum Generated [GeV]")
histo35.GetYaxis().SetTitle("PhiB Simulated [rad]")
histo35.SetTitle("PhiB at Station 1 entry")
histo35.SetStats(0)
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
histo28.GetXaxis().SetTitle("Transverse Momentum Generated [GeV]")
histo28.GetYaxis().SetTitle("PhiB Simulated [rad]")
histo28.SetTitle("PhiB at Station 2 entry")
histo28.SetStats(0)
#fitFunc1=TF1("fitFunc1", "pol1", -100, 100)
#histo28.Fit("fitFunc1", "R", " ")
#fitFunc1.Draw('same')
histo28.DrawCopy('COL')
c28.Print("./plots/"+c28.GetName()+".png")
c28.Update()
input('press enter to exit')


'''



#zad 21
c29 = TCanvas('cPhiBCompSt1', 'PhiB Comb', 600, 600)
c29.SetLeftMargin(0.15)  #Space for printing Y label
histo29 = gROOT.FindObject('hPhiBCompSt1')
histo29.GetYaxis().SetTitle("PhiB Reconstructed [rad]")
histo29.GetXaxis().SetTitle("PhiB Simulated [rad]")
histo29.SetTitle("Comparing PhiB at Station 1 entry")
histo29.SetStats(0)
histo29.DrawCopy('COL')
histo29.GetXaxis().SetRange(-1, 1)



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
histo34.GetYaxis().SetTitle("PhiB Reconstructed [rad]")
histo34.GetXaxis().SetTitle("PhiB Simulated [rad]")
histo34.SetTitle("Comparing PhiB at Station 2 entry")
histo34.SetStats(0)


histo34.DrawCopy('COL')
line = TLine(-1, -1, 1, 1)
line.SetLineColor(2)
#line.SetLineWidth(2)
line.Draw("same")
c34.Print("./plots/"+c34.GetName()+".png")
c34.Update()
input('press enter to exit')



#zad 25.1
c32 = TCanvas('cPhiCompareSt1', 'Phi Comp st 1', 600, 600)
c32.SetLeftMargin(0.15)  #Space for printing Y label
histo32 = gROOT.FindObject('hPhiCompSt1')
histo32.GetYaxis().SetTitle("Phi Reconstructed [rad]")
histo32.GetXaxis().SetTitle("Phi Simulated [rad]")
histo32.SetTitle("Comparing Phi at Station 1 entry")
#not working 
#histo32.SetMarkerStyle(20)
#histo32.SetMarkerSize(10.0) #Point Size
#histo32.Draw("P")
histo32.SetStats(0)
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
histo33.GetYaxis().SetTitle("Phi Reconstructed [rad]")
histo33.GetXaxis().SetTitle("Phi Simulated [rad]")
histo33.SetTitle("Comparing Phi at Station 2 entry")
#histo33.SetMarkerSize(3.0)
histo33.SetStats(0)
histo33.DrawCopy('COL')
line = TLine(-4, -4, 4, 4)
line.SetLineColor(2)
#line.SetLineWidth(2)
line.Draw("same")
c33.Print("./plots/"+c33.GetName()+".png")
c33.Update()
input('press enter to exit')

#zad 28.1
c36 = TCanvas('cDeltaPhiB1', 'Delta PhiB at Station 1 entry', 600, 600)
histo36 = gROOT.FindObject('hDeltaPhiB1')
histo36.SetTitle("Delta PhiB at Station 1 entry")
histo36.GetXaxis().SetTitle("Delta PhiB [rad]")
histo36.GetYaxis().SetTitle("Entries")
#histo36.SetStats(0)
c36.SetLeftMargin(0.15)  #Space for printing Y label
histo36.DrawCopy('COL')
c36.Print("./plots/"+c36.GetName()+".png")
c36.Update()
input('press enter to exit')

#zad 28.2
c37 = TCanvas('cDeltaPhiB2', 'Delta PhiB at Station 2 entry', 600, 600)
histo37 = gROOT.FindObject('hDeltaPhiB2')
histo37.SetTitle("Delta PhiB at Station 2 entry")
histo37.GetXaxis().SetTitle("Delta PhiB [rad]")
histo37.GetYaxis().SetTitle("Entries")
#histo37.SetStats(0)
c37.SetLeftMargin(0.15)  #Space for printing Y label
histo37.DrawCopy('COL')
c37.Print("./plots/"+c37.GetName()+".png")
c37.Update()
input('press enter to exit')

#zad 29.1
c38 = TCanvas('cDeltaPhi1', 'Delta Phi at station 1 entry', 600, 600)
histo38 = gROOT.FindObject('hDeltaPhi1')
histo38.SetTitle("Delta Phi at Station 1 entry")
histo38.GetXaxis().SetTitle("Delta Phi [rad]")
histo38.GetYaxis().SetTitle("Entries")
#histo38.SetStats(0)
c38.SetLeftMargin(0.15)  #Space for printing Y label
histo38.DrawCopy('COL')
c38.Print("./plots/"+c38.GetName()+".png")
c38.Update()
input('press enter to exit')

#zad 29.2
c39 = TCanvas('cDeltaPhi2', 'Delta Phi at station 2 entry', 600, 600)
histo39 = gROOT.FindObject('hDeltaPhi2')
histo39.SetTitle("Delta Phi at Station 2 entry")
histo39.GetXaxis().SetTitle("Delta Phi [rad]")
histo39.GetYaxis().SetTitle("Entries")
#histo39.SetStats(0)
c39.SetLeftMargin(0.15)  #Space for printing Y label
histo39.DrawCopy('COL')
c39.Print("./plots/"+c39.GetName()+".png")
c39.Update()
input('press enter to exit')


#zad 30.1
histo40 = gROOT.FindObject('hDeltaBCodeSt1')
for (code_value) in range(2, 7):
    bin_code = histo40.GetXaxis().FindBin(code_value)
    histo_projection = histo40.ProjectionY(f"histo_projection_{code_value}", bin_code, bin_code)

    c40 = TCanvas(f'cSt1Code{code_value}', f'Delta PhiB at Station 1 code {code_value}', 600, 600)
    histo_projection.GetXaxis().SetTitle(f'DeltaPhiB (Code = {code_value}) [rad]')
    histo_projection.GetYaxis().SetTitle('Entries')
    histo_projection.SetTitle('Delta PhiB at Station 1')
    #histo_projection.SetStats(0)
    histo_projection.Draw("COL")

    c40.Print(f"./plots/{c40.GetName()}.png")
    c40.Update()

print("Finished generating projections and saving plots")

#zad 30.2
histo41 = gROOT.FindObject('hDeltaBCodeSt2')
for (code_value) in range(2, 7):
    bin_code = histo41.GetXaxis().FindBin(code_value)
    histo_projection = histo41.ProjectionY(f"histo_projection_{code_value}", bin_code, bin_code)

    c41 = TCanvas(f'cSt2Code{code_value}', f'Delta PhiB at Station 2 code {code_value}', 600, 600)
    histo_projection.GetXaxis().SetTitle(f'DeltaPhiB (Code = {code_value}) [rad]')
    histo_projection.GetYaxis().SetTitle('Entries')
    histo_projection.SetTitle('Delta PhiB at Station 2')
    #histo_projection.SetStats(0)
    histo_projection.Draw("COL")

    c41.Print(f"./plots/{c41.GetName()}.png")
    c41.Update()

print("Finished generating projections and saving plots")


#zad 31 
c42 = TCanvas('cQuality_Compare', 'Old vs New (HW Base) Quality', 600, 600)
histo42 = gROOT.FindObject('hQuality_Compare')
histo42.SetTitle("Old vs New (HW Base) Quality")
histo42.GetXaxis().SetTitle("Quality Old")
histo42.GetYaxis().SetTitle("Quality New")
histo42.SetStats(0)
c42.SetLeftMargin(0.15)  #Space for printing Y label
histo42.DrawCopy('COL')
c42.Print("./plots/"+c42.GetName()+".png")
c42.Update()
input('press enter to exit')