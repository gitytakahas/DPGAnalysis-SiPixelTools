from officialStyle import officialStyle 
from ROOT import TFile, TTree, TH2F, TCanvas, gROOT, gStyle, TH1F, TLegend
import copy

gROOT.SetBatch(True) 
officialStyle(gStyle) 
gStyle.SetOptTitle(0)

def LegendSettings(leg):
    leg.SetBorderSize(0)
    leg.SetFillColor(10)
    leg.SetLineColor(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetTextFont(42)


#type='zerobias'
#type='random'

layers = [1,2,3,4]
#layers = [1]
xmax = [10, 3, 2, 1]
types = ['']
ladder = [13, 29, 45, 65]
lmax = [6.5, 14.5, 22.5, 32.5]

h_occupancy = {}

for type in types:


    file = TFile('Myroot.root')

    for index in range(1, 2):

        cname = 'canvas_' + str(index)
        canvas = TCanvas(cname, cname, 1200,800)

        canvas.Divide(2,2)
        

        for layer in layers:

            canvas.cd(layer)
            canvas.cd(layer).SetGridx()
            canvas.cd(layer).SetGridy()

            hist = file.Get('hmap_L' + str(layer) + '_id' + str(index))
    
            hist.GetXaxis().SetTitle('Z (mm)')
            hist.GetYaxis().SetTitle('Phi')
            
            hist.Draw('colz')
            
        canvas.SaveAs('plot/clustermap_' + str(index) + '.gif')


