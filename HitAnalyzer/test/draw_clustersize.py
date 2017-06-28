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
    tree = file.Get('cluster_tree')

    h_occupancy_ = []

    for layer in layers:
    
        hname = 'hist_L' + str(layer)
#        hist = TH2F(hname, hname, 56, -28, 28, 10000,0,10000)
        hist = TH2F(hname, hname, 20, -28, 28, 10000,0,10000)

        hist.GetXaxis().SetTitle('Z (mm)')
        hist.GetYaxis().SetTitle('Cluster size')

        tree.Draw("ch:zPos >> " + hname, "subid==1 && layer==" + str(layer))

        cname = 'canvas_' + str(layer)
        canvas = TCanvas(cname)
    
        canvas.SetGridx()
        canvas.SetGridy()
        
        hist.Draw('colz')

    
        hist_occ = hist.ProfileX()
        hist_occ.GetYaxis().SetNdivisions(505)
        hist_occ.Sumw2()
        
        hist_occ.SetLineColor(1)
#        hist_occ.Draw('psame')

        canvas.SaveAs('plot/cluster_L'+str(layer) + '_' + type + '.gif')


        ## zoom

        hname_zoom = 'hist_zoom_L' + str(layer)
        hist_zoom = TH2F(hname_zoom, hname_zoom, 20, -28, 28, 100,0,200)

        hist_zoom.GetXaxis().SetTitle('Z (mm)')
        hist_zoom.GetYaxis().SetTitle('Cluster size')

        tree.Draw("ch:zPos >> " + hname_zoom, "subid==1 && layer==" + str(layer))

        cname_zoom = 'canvas_zoom_' + str(layer)
        canvas_zoom = TCanvas(cname_zoom)
    
        canvas_zoom.SetGridx()
        canvas_zoom.SetGridy()
       
        hist_zoom.Draw('colz')
#        hist_occ.Draw('psame')
        hist.Draw('candlex(10000311) same')
        canvas_zoom.SaveAs('plot/cluster_zoom_L'+str(layer) + '_' + type + '.gif')




#        h_occupancy_.append(copy.deepcopy(hist_zoom))
        h_occupancy_.append(copy.deepcopy(hist))
#        h_occupancy_.append(copy.deepcopy(hist_occ))

    h_occupancy[type] = h_occupancy_


print h_occupancy 
#    LegendSettings(leg,len(hists))

gStyle.SetPadRightMargin(0.1)
gStyle.SetPadLeftMargin(0.18)

types.reverse()

for layer in layers:
    cname = 'occupancy_' + str(layer)
    canvas_layer = TCanvas(cname)

    leg = TLegend(0.5,0.7,0.9,0.9)
    LegendSettings(leg)

    for index, type in enumerate(types):
    
#        h_occupancy[type][layer-1].Scale(1./h_occupancy[type][layer-1].GetSumOfWeights())

        h_occupancy[type][layer-1].SetLineWidth(2)
        h_occupancy[type][layer-1].SetLineColor(index+1)
        h_occupancy[type][layer-1].SetMarkerColor(index+1)
        h_occupancy[type][layer-1].SetLineStyle(index+1)

        h_occupancy[type][layer-1].GetXaxis().SetTitle('Z (mm)')
        h_occupancy[type][layer-1].GetYaxis().SetTitle('Cluster size')
        h_occupancy[type][layer-1].GetYaxis().SetTitleOffset(1.5)
        h_occupancy[type][layer-1].GetYaxis().SetRangeUser(0,200)
        h_occupancy[type][layer-1].SetMaximum(h_occupancy[type][layer-1].GetMaximum()*1.5)
        h_occupancy[type][layer-1].SetMinimum(0)

        if index==0:
            h_occupancy[type][layer-1].Draw('h')
            h_occupancy[type][layer-1].Draw('candlex(10000311)')
#            leg.AddEntry(h_occupancy[type][layer-1], 'Layer'+str(layer), '')
        else:
#            h_occupancy[type][layer-1].Draw('hsame')
            h_occupancy[type][layer-1].Draw('candlex(10000311) same')


        leg.AddEntry(h_occupancy[type][layer-1], type, 'lep')

    leg.Draw()

    canvas_layer.SaveAs('plot/cluster_profile_L' + str(layer) + '.gif')
    
