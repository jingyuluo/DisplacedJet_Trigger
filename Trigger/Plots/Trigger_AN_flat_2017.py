import ROOT
import sys
import numpy
import array
import argparse
from math import sqrt, fabs, sin, exp, log10
from ROOT import TFile, TTree, TChain, TBranch, TH1D, TH1I, TH1F, TH2F, TGraphAsymmErrors, TEfficiency, Math
from ROOT import TLorentzVector
from ROOT.Math import LorentzVector

parser = argparse.ArgumentParser(description="Meausre the trigger efficiency")

parser.add_argument("-f", "--file", default="", help="The path to the output tree of displaced jets analysis")
parser.add_argument("-d", "--data", default=False, action="store_true", help="CMSSW_7_4_X or CMSSW_7_6_X")
parser.add_argument("-o", "--output", default="", help="The name of the output file")
parser.add_argument("-w", "--weight", default="1.0",  help="The reweight factor of the sample")
args = parser.parse_args()

filename = args.file
weight = float(args.weight)
tfile = TFile.Open(filename)
tree = tfile.Get("TriggerNtuple/tree")

outname = args.output
fout = TFile("{0}.root".format(outname), "recreate")

#eff_hlt_HT_disp_430 = TEfficiency("eff_hlt_HT_dis_430", "Trigger Efficiency of HT430_displaceddijet40", 200, 100, 1000)

eff_HT_disp_430_uncorr = TEfficiency("eff_HT_disp_430_uncorr", "Trigger Efficiency of HT430_displaceddijet40", 200, 100, 3000)
eff_HT_disp_430_zoom_uncorr = TEfficiency("eff_HT_dis_430_zoom_uncorr", "Trigger Efficiency of HT430_displaceddijet40", 200, 100, 1000)
eff_HT_disp_430 = TEfficiency("eff_HT_dis_430", "Trigger Efficiency of HT430_displaceddijet40", 200, 100, 3000)
eff_HT_disp_430_zoom = TEfficiency("eff_HT_dis_430_zoom", "Trigger Efficiency of HT430_displaceddijet40", 200, 100, 1000)

xbins=array.array('d', [0, 100, 150, 200, 220, 240, 260, 280, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 820, 840, 860, 880, 900, 920, 940, 960, 980, 1000, 1040, 1080, 1100, 1120, 1140, 1160, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1600, 2000, 3000])
eff_HT_disp_430_variable = TEfficiency("eff_HT_disp_430_variable", "eff_HT_disp_430_variable", 64, xbins)

#eff_HT_disp_430_3p0 = TEfficiency("eff_HT_dis_430_3p0", "Trigger Efficiency of HT430_displaceddijet40_3p0", 200, 100, 1000)

nentries = tree.GetEntries()

passed_evt = 0

for entry in range(nentries):
   
    tree.GetEntry(entry)
    print "Processing Event: ", entry

    HT430_filter = tree.HTfilter_fired

    #recalculate the offline caloHT

    eff_HT_disp_430_uncorr.Fill(bool(HT430_filter), tree.unCorr_CaloHT)
    eff_HT_disp_430_zoom_uncorr.Fill(bool(HT430_filter), tree.unCorr_CaloHT)
    eff_HT_disp_430.Fill(bool(HT430_filter), tree.CaloHT)
    eff_HT_disp_430_variable.Fill(bool(HT430_filter), tree.CaloHT)
    eff_HT_disp_430_zoom.Fill(bool(HT430_filter), tree.CaloHT)

gra_HTdisp_430 = eff_HT_disp_430.CreateGraph()
print passed_evt
#fout.WriteTObject(gra_HTdisp_430, "gra_HTdisp_430")
fout.Write()
fout.Close()

