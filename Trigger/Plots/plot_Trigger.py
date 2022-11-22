import sys, os
import argparse

sys.path.append('/Users/jingyuluo/Downloads/root_dir/lib')

import ROOT

from ROOT import TFile, TCanvas, TH1F
from ROOT import gStyle
from ROOT import kGreen, kRed, kBlue

gStyle.SetOptStat(ROOT.kFALSE)
parser=argparse.ArgumentParser()

parser.add_argument("-f", "--file", help="The ROOT file for the histogram")
parser.add_argument("-l", "--label", help="The label for the output file")

args = parser.parse_args()

filename=args.file

can = TCanvas("can", "can", 800, 800)
can.cd()
can.SetTickx()
can.SetTicky()
can.SetGrid(1,1)

tfile = TFile(filename)
eff_HT430 = tfile.Get("eff_HT_disp_430_variable")

gra_HT430 = eff_HT430.CreateGraph()
gra_HT430.SetTitle("")
gra_HT430.SetMarkerColor(kBlue)
gra_HT430.SetMarkerStyle(ROOT.kFullCircle)
gra_HT430.SetMarkerSize(0.75)
gra_HT430.GetXaxis().SetTitle("Offline Calo H_{T} [GeV]")
gra_HT430.GetXaxis().SetRangeUser(0,3000)
gra_HT430.GetYaxis().SetTitle("Calo H_{T} Filter Efficiency")
gra_HT430.Draw("APE")
text=ROOT.TLatex(0.65, 0.92, "2018 (13TeV)")
text.SetNDC()
text.SetTextFont(62)
text.SetTextSize(0.05)
text2=ROOT.TLatex(0.15, 0.92, "CMS #bf{#scale[0.75]{#it{Preliminary}}}")
text2.SetNDC()
text2.SetTextSize(0.05)
text2.SetTextFont(62)
text.Draw("SAME")
text2.Draw("SAME")
leg=ROOT.TLegend(0.6,0.45, 0.85,0.6)
leg.SetLineColor(0)
leg.SetFillColorAlpha(0,0)

leg.AddEntry(gra_HT430, "#bf{Run2018}", "pe")
leg.Draw("SAME")

can.SaveAs("Run2018All_hltCaloHT_Eff.png")
can.SaveAs("Run2018All_hltCaloHT_Eff.pdf")
can.SaveAs("Run2018All_hltCaloHT_Eff.eps")
can.SaveAs("Run2018All_hltCaloHT_Eff.C")
