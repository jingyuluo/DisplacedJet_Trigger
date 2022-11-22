import sys, os
import argparse

sys.path.append('/Users/jingyuluo/Downloads/root_dir/lib')

import ROOT

from ROOT import TFile, TCanvas, TH1F
from ROOT import gStyle
from ROOT import kGreen, kRed, kBlue, kBlack

gStyle.SetOptStat(ROOT.kFALSE)
parser=argparse.ArgumentParser()

parser.add_argument("-d1", "--data1", default="test_Run2022C_SM_trigger.root", help="The ROOT file for the 2022 trigger performance")
parser.add_argument("-d2", "--data2", default="SingleMuon_Run2018_merged.root", help="The ROOT file for the 2018 trigger performance")
parser.add_argument("-m", "--mc", default="SingleMuon_QCD_merged.root",help="The ROOT file for the MC trigger performance")

args = parser.parse_args()


filename1 = args.data1
filename2 = args.data2
filemc = args.mc

can = TCanvas("can", "can", 800, 800)
can.cd()
can.SetTickx()
can.SetTicky()
can.SetGrid(1,1)

tfile1 = TFile(filename1)
tfile2 = TFile(filename2)
tfilemc = TFile(filemc)
eff_HT430 = tfile1.Get("eff_HT_dis_430_zoom")
eff_HT430_2018  = tfile2.Get("eff_HT_dis_430_zoom")
eff_HT430_mc = tfilemc.Get("eff_HT_dis_430_zoom")

gra_HT430 = eff_HT430.CreateGraph()
gra_HT430_2018 = eff_HT430_2018.CreateGraph()
gra_HT430_mc = eff_HT430_mc.CreateGraph()

gra_HT430.SetTitle("")
gra_HT430.SetMarkerColor(kBlack)
gra_HT430.SetMarkerStyle(ROOT.kFullCircle)
gra_HT430.SetMarkerSize(0.75)
gra_HT430.GetXaxis().SetTitle("Offline Calo H_{T} [GeV]")
gra_HT430.GetXaxis().SetRangeUser(0,3000)
gra_HT430.GetYaxis().SetTitle("Calo H_{T} Filter Efficiency")

gra_HT430_2018.SetMarkerColor(kBlue)
gra_HT430_2018.SetMarkerStyle(ROOT.kFullCircle)
gra_HT430_2018.SetMarkerSize(0.75)
gra_HT430_2018.SetLineColor(kBlue)

gra_HT430_mc.SetMarkerColor(kRed)
gra_HT430_mc.SetMarkerStyle(ROOT.kFullCircle)
gra_HT430_mc.SetMarkerSize(0.75)
gra_HT430_mc.SetLineColor(kRed)

gra_HT430.Draw("APE")
gra_HT430_2018.Draw("PE SAME")
gra_HT430_mc.Draw("PE SAME")

text=ROOT.TLatex(0.65, 0.92, "2022")
text.SetNDC()
text.SetTextFont(62)
text.SetTextSize(0.05)
text2=ROOT.TLatex(0.15, 0.92, "CMS #bf{#scale[0.75]{#it{Preliminary}}}")
text2.SetNDC()
text2.SetTextSize(0.05)
text2.SetTextFont(62)
text.Draw("SAME")
text2.Draw("SAME")
leg=ROOT.TLegend(0.52,0.45, 0.88,0.6)
leg.SetLineColor(0)
leg.SetFillColorAlpha(0,0)

leg.AddEntry(gra_HT430, "Run2022", "pe")
leg.AddEntry(gra_HT430_2018, "Run2018", "pe")
leg.AddEntry(gra_HT430_mc, "Run2 QCD MC", "pe")

leg.Draw("SAME")

can.SaveAs("Run2022C_hltCaloHT_Eff.png")
can.SaveAs("Run2022C_hltCaloHT_Eff.pdf")
can.SaveAs("Run2022C_hltCaloHT_Eff.eps")
can.SaveAs("Run2022C_hltCaloHT_Eff.C")
