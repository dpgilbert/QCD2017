import pyRootPlotMaker as ppm
import ROOT
import sys

ROOT.gROOT.SetBatch(True)

# if len(sys.argv) > 1 and sys.argv[1] == "mode":

name16="~/MT2AnalysisEdit/MyCode/QCDcheck_16vs17/output/data_2016.root"
name17="~/2017/MT2AnalysisSep13/QCDcheck/output/data_2017.root"

f16 = ROOT.TFile(name16)
f17 = ROOT.TFile(name17)

hists=["h_mt2","h_ht","h_ht_mt2cut","h_j1pt","h_j2pt","h_met_pt","h_met_phi","h_jet1_phi", "h_mt2_hi","h_ht_hi","h_ht_mt2cut_hi","h_j1pt_hi","h_j2pt_hi","h_met_pt_hi","h_met_phi_hi","h_jet1_phi_hi"]

for hist in hists:

    h16=f16.Get(hist)
    h17=f17.Get(hist)

    h16.Scale(11.379/35.9)

    h_16_vec = [h16]
    h_16_names = ["Data (2016, Scaled)"]

    if hist.find("_hi") == -1:
        htmet = "MET > 30 GeV; If H_{T} < 1200, MET > 300 GeV"
    else:
        htmet = "MET > 300 GeV"

    y=[0,2]
    if hist.find("h_mt2") != -1:
        desc="M_{T2} (H_{T} > 300 GeV)"
        xtitle="M_{T2}"
        xunit="GeV"
        log=True
        y=[0.5,1.5]
    elif hist.find("h_ht") != -1:
        desc="H_{T} (Inclusive in M_{T2})"
        xtitle="H_{T}"
        xunit="GeV"
        log=True
        y=[0,2]
    elif hist.find("h_ht_mt2cut") != -1:
        desc="H_{T} (M_{T2} > 200 GeV)"
        xtitle="H_{T}"
        xunit="GeV"
        log=True
        y=[0,3]
    elif hist.find("h_j1pt") != -1:
        desc="p_{T}^{Jet 1} (H_{T} > 300 GeV)"
        xtitle=desc
        xunit="GeV"
        log=True
        y=[0.5,3]
    elif hist.find("h_j2pt") != -1:
        desc="p_{T}^{Jet 2} (H_{T} > 300 GeV)"
        xtitle=desc
        xunit="GeV"
        log=True
        y=[0,2]
    elif hist.find("h_met_pt") != -1:
        desc="MET (H_{T} > 1200 GeV)"
        xtitle="MET"
        xunit="GeV"
        log=True
        y=[0.5,5]
    elif hist.find("h_met_phi") != -1:
        desc="#phi_{MET} (H_{T} > 1200 GeV)"
        xtitle="#phi_{MET}"
        xunit="rad"
        log=False
        y=[0.5,3]
    elif hist.find("h_jet1_phi") != -1:
        desc="#phi_{jet 1} (H_{T} > 1200 GeV)"
        xtitle="#phi_{Jet 1}"
        xunit="rad"
        log=False
        y=[0.5,3]

    jet = "All-hadronic, 2+ Jets"

    metmht = "|MET-MHT|/MET < 0.5"

    phimin = "QCD Enriched (#Delta#phi(Jet_{1...4},MET)_{min} < 0.3)"

    subtitles = [desc,jet,htmet,metmht,phimin]

    outputFiles = ["PNGs/"+hist+".png","PDFs/"+hist+".pdf"]

    for outputFile in outputFiles:
        ppm.plotDataMC(h_16_vec,h_16_names,h_data=h17,dataTitle="Data (2017)",isLog = log,ratioTitle="2017 / 2016", subtitles=subtitles, subtitleSize=0.02,
                           xAxisTitle=xtitle, xAxisUnit=xunit, yRangeUserRatio=y, markerSize=1.1,lumi=11.379,lumiUnit="fb",energy=13,saveAs=outputFile)

