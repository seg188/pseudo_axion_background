from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel
import ROOT as root

file = root.TFile.Open("~/hex/sp_data.root")
data = file.Get("sp_data_MassPi0_PT_0_side");
data.Scale(1/data.Integral());

datax = []
datay = []
for k in range(data.GetNbinsX()):
	datax.append([])
	datax[k].append(data.GetBinCenter(k));
	datay.append(data.GetBinContent(k));

kernel = DotProduct() + WhiteKernel()
gpr = GaussianProcessRegressor().fit(datax, datay)
gpr.score(datax, datay)

ndata = 100;
minv = 0.0
maxv = 5.0

graph = root.TGraph(ndata)
_sigma_1_graph = root.TGraph(2*ndata)

 

for x in range(ndata):
	xv = minv + x*(maxv-minv)/ndata
	xv2 = maxv - minv + (x + 1)*(maxv-minv)/ndata

	y, dy =  gpr.predict([[xv]], return_std=True)
	graph.SetPoint(x, xv, y)
	_sigma_1_graph.SetPoint(x, xv, y + dy)
	_sigma_1_graph.SetPoint(x + ndata, xv2, y-dy)


_sigma_1_graph.SetFillStyle(3005);
_sigma_1_graph.SetFillColor(2);

canv = root.TCanvas()
data.Draw()
graph.Draw("SAME")
#_sigma_1_graph.Draw("SAME")

canv.Print("test.png", ".png")