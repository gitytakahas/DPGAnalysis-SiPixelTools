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
types = ['random', 'zerobias']
ladder = [13, 29, 45, 65]
lmax = [6.5, 14.5, 22.5, 32.5]

h_occupancy = {}

for type in types:


    file = TFile('Myroot_' + type + '.root')
    tree = file.Get('pixel_tree')

    hnorm = file.Get('total')
    norm = hnorm.GetBinContent(1)

    h_occupancy_ = []

    for layer in layers:
    
        hname = 'hist_L' + str(layer)
        hist = TH2F(hname, hname, 72, -4.5, 4.5, 2*ladder[layer-1], -lmax[layer-1], lmax[layer-1])

        hist.GetXaxis().SetTitle('module')
        hist.GetYaxis().SetTitle('ladder')

        #    tree.Draw("tmp2:tmp1 >> " + hname, "(subid==1 && layer==" + str(layer) + ")*ch")
        tree.Draw("tmp2:tmp1 >> " + hname, "(subid==1 && layer==" + str(layer) + ")")

        cname = 'canvas_' + str(layer)
        canvas = TCanvas(cname)
    
        canvas.SetGridx()
        canvas.SetGridy()
        
        hist.Scale(1./norm)
        hist.Draw('colz')
        canvas.SaveAs('plot/L'+str(layer) + '_' + type + '.gif')


        oname = 'hist_occ_L' + str(layer)

        hist_occ = TH1F(oname, oname, 50,0,xmax[layer-1])
        hist_occ.GetYaxis().SetNdivisions(505)
        hist_occ.Sumw2()

        for ix in range(1, hist.GetXaxis().GetNbins()+1):
            for iy in range(1, hist.GetYaxis().GetNbins()+1):
                val = hist.GetBinContent(ix, iy)
                print ix, iy, val
                
                if val!=0:
                    hist_occ.Fill(val)



        h_occupancy_.append(copy.deepcopy(hist_occ))

    h_occupancy[type] = h_occupancy_


print h_occupancy 
#    LegendSettings(leg,len(hists))

gStyle.SetPadRightMargin(0.1)
gStyle.SetPadLeftMargin(0.18)


for layer in layers:
    cname = 'occupancy_' + str(layer)
    canvas_layer = TCanvas(cname)

    leg = TLegend(0.5,0.7,0.9,0.9)
    LegendSettings(leg)

    for index, type in enumerate(types):
    
        h_occupancy[type][layer-1].Scale(1./h_occupancy[type][layer-1].GetSumOfWeights())

        h_occupancy[type][layer-1].SetLineWidth(2)
        h_occupancy[type][layer-1].SetLineColor(index+1)
        h_occupancy[type][layer-1].SetMarkerColor(index+1)
        h_occupancy[type][layer-1].SetLineStyle(index+1)

        h_occupancy[type][layer-1].GetXaxis().SetTitle('Pixel hit occupancy / module / event')
        h_occupancy[type][layer-1].GetYaxis().SetTitle('Normalized')
        h_occupancy[type][layer-1].GetYaxis().SetTitleOffset(1.5)
        h_occupancy[type][layer-1].SetMinimum(0)

        if index==0:
            h_occupancy[type][layer-1].Draw('h')
            leg.AddEntry(h_occupancy[type][layer-1], 'Layer'+str(layer), '')
        else:
            h_occupancy[type][layer-1].Draw('hsame')


        legname = type + ' (' + '{0:.2f}'.format(h_occupancy[type][layer-1].GetMean())  + ', {0:.2f}'.format(h_occupancy[type][layer-1].GetRMS()) +')'

        leg.AddEntry(h_occupancy[type][layer-1], legname, 'lep')

    leg.Draw()

    canvas_layer.SaveAs('plot/Occupancy_L' + str(layer) + '.gif')
    
