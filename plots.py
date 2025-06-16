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



# #przyklad
# c0 = TCanvas('cLicExample','cLicExample',600,600)
# histo0 = gROOT.FindObject('hLicExample')
# histo0.DrawCopy('box text')
# c0.Print("./plots/"+c0.GetName()+".png")
# c0.Update()
# input('press enter to exit')


# #zad 27
# c35 = TCanvas('cPhiB_st1', 'PhiB(Pt) St1', 600, 600)
# c35.SetLeftMargin(0.15)  #Space for printing Y label
# c35.SetLogy(True)
# histo35 = gROOT.FindObject('hPhiB_st1')
# histo35.GetXaxis().SetTitle("Transverse Momentum Generated [GeV]")
# histo35.GetYaxis().SetTitle("PhiB Simulated [rad]")
# histo35.SetTitle("PhiB at Station 1 entry")
# histo35.SetStats(0)
# #fitFunc1=TF1("fitFunc1", "pol1", -100, 100)
# #histo28.Fit("fitFunc1", "R", " ")
# #fitFunc1.Draw('same')
# histo35.DrawCopy('COL')
# c35.Print("./plots/"+c35.GetName()+".png")
# c35.Update()
# input('press enter to exit')


# #zad 20
# c28 = TCanvas('cPhiB_st2', 'PhiB(Pt) st 2', 600, 600)
# c28.SetLeftMargin(0.15)  #Space for printing Y label
# c28.SetLogy(True)
# histo28 = gROOT.FindObject('hPhiB_st2')
# histo28.GetXaxis().SetTitle("Transverse Momentum Generated [GeV]")
# histo28.GetYaxis().SetTitle("PhiB Simulated [rad]")
# histo28.SetTitle("PhiB at Station 2 entry")
# histo28.SetStats(0)
# #fitFunc1=TF1("fitFunc1", "pol1", -100, 100)
# #histo28.Fit("fitFunc1", "R", " ")
# #fitFunc1.Draw('same')
# histo28.DrawCopy('COL')
# c28.Print("./plots/"+c28.GetName()+".png")
# c28.Update()
# input('press enter to exit')






# #zad 21
# c29 = TCanvas('cPhiBCompSt1', 'PhiB Comb', 600, 600)
# c29.SetLeftMargin(0.15)  #Space for printing Y label
# histo29 = gROOT.FindObject('hPhiBCompSt1')
# histo29.GetYaxis().SetTitle("PhiB Reconstructed [rad]")
# histo29.GetXaxis().SetTitle("PhiB Simulated [rad]")
# histo29.SetTitle("Comparing PhiB at Station 1")
# histo29.SetStats(0)
# histo29.DrawCopy('COL')
# histo29.GetXaxis().SetRange(-1, 1)



# line = TLine(-1, -1, 1, 1)
# line.SetLineColor(2)
# #line.SetLineWidth(2)
# line.Draw("same")

# c29.Print("./plots/"+c29.GetName()+".png")
# c29.Update()
# input('press enter to exit')



# #zad 26
# c34 = TCanvas('cPhiBCompSt2', 'PhiB Comb', 600, 600)
# c34.SetLeftMargin(0.15)  #Space for printing Y label
# histo34 = gROOT.FindObject('hPhiBCompSt2')
# histo34.GetYaxis().SetTitle("PhiB Reconstructed [rad]")
# histo34.GetXaxis().SetTitle("PhiB Simulated [rad]")
# histo34.SetTitle("Comparing PhiB at Station 2")
# histo34.SetStats(0)


# histo34.DrawCopy('COL')
# line = TLine(-1, -1, 1, 1)
# line.SetLineColor(2)
# #line.SetLineWidth(2)
# line.Draw("same")
# c34.Print("./plots/"+c34.GetName()+".png")
# c34.Update()
# input('press enter to exit')



# #zad 25.1
# c32 = TCanvas('cPhiCompareSt1', 'Phi Comp st 1', 600, 600)
# c32.SetLeftMargin(0.15)  #Space for printing Y label
# histo32 = gROOT.FindObject('hPhiCompSt1')
# histo32.GetYaxis().SetTitle("Phi Reconstructed [rad]")
# histo32.GetXaxis().SetTitle("Phi Simulated [rad]")
# histo32.SetTitle("Comparing Phi at Station 1")
# # not working 
# # histo32.SetMarkerStyle(20)
# # histo32.SetMarkerSize(10.0) #Point Size
# # histo32.Draw("P")
# histo32.SetStats(0)
# histo32.DrawCopy('COL')
# line = TLine(-4, -4, 4, 4)
# line.SetLineColor(2)
# #line.SetLineWidth(2)
# line.Draw("same")
# c32.Print("./plots/"+c32.GetName()+".png")
# c32.Update()
# input('press enter to exit')

# #zad 25.2
# c33 = TCanvas('cPhiCompareSt2', 'Phi Comp st 2', 600, 600)
# c33.SetLeftMargin(0.15)  #Space for printing Y label
# histo33 = gROOT.FindObject('hPhiCompSt2')
# histo33.GetYaxis().SetTitle("Phi Reconstructed [rad]")
# histo33.GetXaxis().SetTitle("Phi Simulated [rad]")
# histo33.SetTitle("Comparing Phi at Station 2")
# #histo33.SetMarkerSize(3.0)
# histo33.SetStats(0)
# histo33.DrawCopy('COL')
# line = TLine(-4, -4, 4, 4)
# line.SetLineColor(2)
# #line.SetLineWidth(2)
# line.Draw("same")
# c33.Print("./plots/"+c33.GetName()+".png")
# c33.Update()
# input('press enter to exit')



#zad 28.1
c36 = TCanvas('cDeltaPhiB1', 'Delta PhiB at Station 1 entry', 600, 600)
histo36 = gROOT.FindObject('hDeltaPhiB1')
histo36.SetTitle("Delta PhiB at Station 1")
histo36.GetXaxis().SetTitle("Delta PhiB [rad]")
histo36.GetYaxis().SetTitle("Entries")
#histo36.SetStats(0)
c36.SetLeftMargin(0.15)  #Space for printing Y label
histo36.DrawCopy('COL')
c36.Print("./plots/"+c36.GetName()+".png")
c36.Update()
input('press enter to exit')

# #zad 28.2
# c37 = TCanvas('cDeltaPhiB2', 'Delta PhiB at Station 2 entry', 600, 600)
# histo37 = gROOT.FindObject('hDeltaPhiB2')
# histo37.SetTitle("Delta PhiB at Station 2")
# histo37.GetXaxis().SetTitle("Delta PhiB [rad]")
# histo37.GetYaxis().SetTitle("Entries")
# #histo37.SetStats(0)
# c37.SetLeftMargin(0.15)  #Space for printing Y label
# histo37.DrawCopy('COL')
# c37.Print("./plots/"+c37.GetName()+".png")
# c37.Update()
# input('press enter to exit')

#zad 29.1
#calculate std dev. 
c38 = TCanvas('cDeltaPhi1', 'Delta Phi at station 1 entry', 600, 600)
histo38 = gROOT.FindObject('hDeltaPhi1')
histo38.SetTitle("Delta Phi at Station 1")
# c38.SetLogx(True)
histo38.GetXaxis().SetTitle("Delta Phi [rad]")
histo38.GetYaxis().SetTitle("Entries")
#histo38.SetStats(0)
histo38.GetXaxis().SetRangeUser(-0.01, 0.01)
histo38.GetXaxis().SetNdivisions(507) #7 labels, 5 subdivisions
c38.SetLeftMargin(0.15)  #Space for printing Y label
histo38.DrawCopy('COL')
c38.Print("./plots/"+c38.GetName()+".png")
c38.Update()
input('press enter to exit')

# # # Manual calculation of mean
# # sum_x = 0
# # sum_w = 0
# # nbins = histo38.GetNbinsX()
# # # for i in range(1, nbins + 1):  # bin numbering starts at 1
# # for i in range(951, 1051):  # bin numbering starts at 1
# #     x = histo38.GetBinCenter(i)
# #     w = histo38.GetBinContent(i)
# #     sum_x += x * w
# #     sum_w += w
# # mean = sum_x / sum_w
# # # Now manual calculation of std dev
# # sum_squared_diff = 0
# # # for i in range(1, nbins + 1):
# # for i in range(951, 1051):  # bin numbering starts at 1
# #     x = histo38.GetBinCenter(i)
# #     w = histo38.GetBinContent(i)
# #     sum_squared_diff += w * (x - mean) ** 2
# # std_dev = math.sqrt(sum_squared_diff / sum_w)
# # print(f"Manual mean = {mean}")
# # print(f"Manual standard deviation = {std_dev}")

# #zad 29.2
# c39 = TCanvas('cDeltaPhi2', 'Delta Phi at station 2 entry', 600, 600)
# histo39 = gROOT.FindObject('hDeltaPhi2')
# histo39.SetTitle("Delta Phi at Station 2")
# histo39.GetXaxis().SetTitle("Delta Phi [rad]")
# histo39.GetYaxis().SetTitle("Entries")
# #histo39.SetStats(0)
# histo39.GetXaxis().SetRangeUser(-0.01, 0.01)
# histo39.GetXaxis().SetNdivisions(507) #7 labels, 5 subdivisions
# c39.SetLeftMargin(0.15)  #Space for printing Y label
# histo39.DrawCopy('COL')
# c39.Print("./plots/"+c39.GetName()+".png")
# c39.Update()
# input('press enter to exit')

# # # Manual calculation of mean
# # sum_x = 0
# # sum_w = 0
# # nbins = histo39.GetNbinsX()
# # # for i in range(1, nbins + 1):  # bin numbering starts at 1
# # for i in range(951, 1051):  # bin numbering starts at 1 (-0.0005; 0.0005) -> Range(951; 1051) for 2000 bins
# #     x = histo39.GetBinCenter(i)
# #     w = histo39.GetBinContent(i)
# #     sum_x += x * w
# #     sum_w += w
# # mean = sum_x / sum_w
# # # Now manual calculation of std dev
# # sum_squared_diff = 0
# # # for i in range(1, nbins + 1):
# # for i in range(951, 1051):  # bin numbering starts at 1
# #     x = histo39.GetBinCenter(i)
# #     w = histo39.GetBinContent(i)
# #     sum_squared_diff += w * (x - mean) ** 2
# # std_dev = math.sqrt(sum_squared_diff / sum_w)
# # print(f"Manual mean = {mean}")
# # print(f"Manual standard deviation = {std_dev}")


# #zad 30.1
# histo40 = gROOT.FindObject('hDeltaBCodeSt1')
# for (code_value) in range(2, 7):
#     bin_code = histo40.GetXaxis().FindBin(code_value)
#     histo_projection = histo40.ProjectionY(f"histo_projection_{code_value}", bin_code, bin_code)

#     c40 = TCanvas(f'cSt1Code{code_value}', f'Delta PhiB at Station 1 code {code_value}', 600, 600)
#     c40.SetLeftMargin(0.15)  #Space for printing Y label
#     histo_projection.GetXaxis().SetTitle(f'DeltaPhiB (Code = {code_value}) [rad]')
#     histo_projection.GetYaxis().SetTitle('Entries')
#     histo_projection.SetTitle('Delta PhiB at Station 1')


#     if code_value == 2:
#         histo_projection.Rebin(8)
#     if code_value == 3:
#         histo_projection.Rebin(7)
#     if code_value == 4:
#         histo_projection.Rebin(4)
#         histo_projection.GetXaxis().SetRangeUser(-0.1, 0.1)
#     if code_value == 5:
#         histo_projection.GetXaxis().SetRangeUser(-0.05, 0.05)
#     if code_value == 6:
#         histo_projection.GetXaxis().SetRangeUser(-0.02, 0.02)
#         histo_projection.SetNdivisions(507) #7 labels, 5 subdivisions

#     #histo_projection.SetStats(0)
#     histo_projection.Draw("COL")

#     c40.Print(f"./plots/{c40.GetName()}.png")
#     c40.Update()

# print("Finished generating projections and saving plots")

# #zad 30.2
# histo41 = gROOT.FindObject('hDeltaBCodeSt2')
# for (code_value) in range(2, 7):
#     bin_code = histo41.GetXaxis().FindBin(code_value)
#     histo_projection = histo41.ProjectionY(f"histo_projection_{code_value}", bin_code, bin_code)

#     c41 = TCanvas(f'cSt2Code{code_value}', f'Delta PhiB at Station 2 code {code_value}', 600, 600)
#     c41.SetLeftMargin(0.15)  #Space for printing Y label
#     histo_projection.GetXaxis().SetTitle(f'DeltaPhiB (Code = {code_value}) [rad]')
#     histo_projection.GetYaxis().SetTitle('Entries')
#     histo_projection.SetTitle('Delta PhiB at Station 2')

#     if code_value == 2:
#         histo_projection.Rebin(8)
#     if code_value == 3:
#         histo_projection.Rebin(7)
#     if code_value == 4:
#         histo_projection.Rebin(4)
#         histo_projection.GetXaxis().SetRangeUser(-0.1, 0.1)
#     if code_value == 5:
#         histo_projection.GetXaxis().SetRangeUser(-0.05, 0.05)
#     if code_value == 6:
#         histo_projection.GetXaxis().SetRangeUser(-0.02, 0.02)
#         histo_projection.SetNdivisions(507) #7 labels, 5 subdivisions

#     #histo_projection.SetStats(0)
#     histo_projection.Draw("COL")

#     c41.Print(f"./plots/{c41.GetName()}.png")
#     c41.Update()

# print("Finished generating projections and saving plots")


# #zad 31 
# c42 = TCanvas('cQuality_Compare', 'Old vs New (HW Base) Quality', 600, 600)
# histo42 = gROOT.FindObject('hQuality_Compare')
# histo42.SetTitle("Quality Codes Comparison in Phase-1 Convention")
# histo42.GetXaxis().SetTitle("Quality Codes in Phase-1")
# histo42.GetYaxis().SetTitle("Quality Codes in Phase-2")

# # histo42.GetXaxis().SetRangeUser(0, 9)
# # histo42.GetYaxis().SetRangeUser(-1, 8)

# histo42.SetStats(0)
# c42.SetLeftMargin(0.15)  #Space for printing Y label
# histo42.DrawCopy('COL')
# histo42.DrawCopy('text, same')
# c42.Print("./plots/"+c42.GetName()+".png")
# c42.Update()
# input('press enter to exit')

# #zad 32.1
# histo43 = gROOT.FindObject('hDeltaCodeSt1')
# for (code_value) in range(2, 7):
#     bin_code = histo43.GetXaxis().FindBin(code_value)
#     histo_projection = histo43.ProjectionY(f"histo_projection_{code_value}", bin_code, bin_code)

#     c43 = TCanvas(f'cPhiSt1Code{code_value}', f'Delta Phi at Station 1 code {code_value}', 600, 600)
#     c43.SetLeftMargin(0.15)  #Space for printing Y label
#     histo_projection.GetXaxis().SetTitle(f'DeltaPhi (Code = {code_value}) [rad]')
#     histo_projection.GetYaxis().SetTitle('Entries')
#     histo_projection.SetTitle('Delta Phi at Station 1')

#     if code_value == 2:
#         histo_projection.Rebin(4)
#         histo_projection.GetXaxis().SetRangeUser(-0.04, 0.04)
#     if code_value == 3:
#         histo_projection.Rebin(8)
#         histo_projection.GetXaxis().SetRangeUser(-0.02, 0.02)
#     if code_value == 4:
#         histo_projection.Rebin()
#         histo_projection.GetXaxis().SetRangeUser(-0.01, 0.01)
#         histo_projection.SetNdivisions(507) #7 labels, 5 subdivisions
#     if code_value == 5:
#         histo_projection.GetXaxis().SetRangeUser(-0.005, 0.005)
#         histo_projection.SetNdivisions(507) #7 labels, 5 subdivisions
#     if code_value == 6:
#         histo_projection.GetXaxis().SetRangeUser(-0.002, 0.002)
#         histo_projection.SetNdivisions(507) #7 labels, 5 subdivisions


#     #histo_projection.SetStats(0)
#     histo_projection.Draw("COL")

#     c43.Print(f"./plots/{c43.GetName()}.png")
#     c43.Update()

# print("Finished generating projections and saving plots")

# #zad 32.2
# histo44 = gROOT.FindObject('hDeltaCodeSt2')
# for (code_value) in range(2, 7):
#     bin_code = histo44.GetXaxis().FindBin(code_value)
#     histo_projection = histo44.ProjectionY(f"histo_projection_{code_value}", bin_code, bin_code)

#     c44 = TCanvas(f'cPhiSt2Code{code_value}', f'Delta Phi at Station 2 code {code_value}', 600, 600)
#     c44.SetLeftMargin(0.15)  #Space for printing Y label
#     histo_projection.GetXaxis().SetTitle(f'DeltaPhi (Code = {code_value}) [rad]')
#     histo_projection.GetYaxis().SetTitle('Entries')
#     histo_projection.SetTitle('Delta Phi at Station 2')

#     if code_value == 2:
#         histo_projection.Rebin(4)
#         histo_projection.GetXaxis().SetRangeUser(-0.04, 0.04)
#     if code_value == 3:
#         histo_projection.Rebin(8)
#         histo_projection.GetXaxis().SetRangeUser(-0.02, 0.02)
#     if code_value == 4:
#         histo_projection.Rebin()
#         histo_projection.GetXaxis().SetRangeUser(-0.01, 0.01)
#         histo_projection.SetNdivisions(507) #7 labels, 5 subdivisions
#     if code_value == 5:
#         histo_projection.GetXaxis().SetRangeUser(-0.005, 0.005)
#         histo_projection.SetNdivisions(507) #7 labels, 5 subdivisions
#     if code_value == 6:
#         histo_projection.GetXaxis().SetRangeUser(-0.002, 0.002)
#         histo_projection.SetNdivisions(507) #7 labels, 5 subdivisions

#     #histo_projection.SetStats(0)
#     histo_projection.Draw("COL")

#     c44.Print(f"./plots/{c44.GetName()}.png")
#     c44.Update()

# print("Finished generating projections and saving plots")



# c48 = TCanvas('cQualityInEvent', 'Code distribution', 600, 600)
# c48.SetLeftMargin(0.15)  #Space for printing Y label
# histo48 = gROOT.FindObject('hQualityInEvent')
# histo48.GetXaxis().SetTitle("Quality Code in Phase-1 convention")
# histo48.GetYaxis().SetTitle("Entries")
# histo48.GetYaxis().SetRangeUser(0, 45000)
# histo48.SetTitle("Quality Code Distribution in Phase-2")
# histo48.SetStats(0)
# histo48.SetFillColor(9)
# histo48.SetNdivisions(12) #7 labels, 5 subdivisions 12+2*100 = 212
# # histo48.DrawCopy("PFC TEXT")
# histo48.DrawCopy("HIST TEXT")
# c48.Print("./plots/"+c48.GetName()+".png")
# c48.Update()
# input('press enter to exit')

# c49 = TCanvas('cQualityInEvent_Leg', 'Code distribution', 600, 600)
# c49.SetLeftMargin(0.15)  #Space for printing Y label
# histo49 = gROOT.FindObject('hQualityInEvent_Leg')
# histo49.GetXaxis().SetTitle("Quality Code in Phase-1 convention")
# histo49.GetYaxis().SetTitle("Entries")
# histo49.GetYaxis().SetRangeUser(0, 45000)
# histo49.SetTitle("Quality Code Distribution in Phase-1")
# histo49.SetStats(0)
# histo49.SetNdivisions(12) #7 labels, 5 subdivisions 12+2*100 = 212
# histo49.SetFillColor(9)
# histo49.DrawCopy("HIST TEXT")
# c49.Print("./plots/"+c49.GetName()+".png")
# c49.Update()
# input('press enter to exit')
