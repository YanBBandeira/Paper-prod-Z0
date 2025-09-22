import ROOT

# 1. Read data from the .dat file
data_points = []
try:
    with open("../Output/dsig_dm.dat", "r") as f:
        for line in f:
            if line.strip() and not line.strip().startswith("#"): # Skip empty lines and comments
                x, y = map(float, line.strip().split())
                data_points.append((x, y))
except FileNotFoundError:
    print("Error: example.dat not found.")
    exit(1)


# 2. Prepare data for TGraph
x_values = ROOT.std.vector('double')()
y_values = ROOT.std.vector('double')()

for x, y in data_points:
    x_values.push_back(x)
    y_values.push_back(y)

# 3. Create and configure TGraph
graph = ROOT.TGraph(len(x_values), x_values.data(), y_values.data())
#construct the TGraph entry (numebr of points, x array, y array)
graph.SetTitle("Data from .dat file;X-axis;Y-axis")
graph.SetMarkerStyle(20) # Circle marker
graph.SetMarkerColor(ROOT.kBlue)
graph.SetLineColor(ROOT.kBlue)
graph.SetLineWidth(2)
graph.SetMarkerSize(1)
graph.GetXaxis().SetTitleSize(0.05)
graph.GetYaxis().SetTitleSize(0.05)
graph.GetXaxis().SetLabelSize(0.04) 
graph.GetYaxis().SetLabelSize(0.04)
graph.GetXaxis().SetTitleOffset(1.2)
graph.GetYaxis().SetTitleOffset(1.4)
graph.GetYaxis().SetRangeUser(0, max(y_values)*1.2) # Set Y-axis range
graph.GetXaxis().SetRangeUser(min(x_values), max(x_values)) # Set X-axis range
graph.SetXaxisTitle("Invariant Mass (GeV)")
graph.SetYaxisTitle("#frac{d#sigma}{dm} (pb/GeV)")

# 4. Draw and display
c = ROOT.TCanvas("c", "Plot from .dat", 800, 600)
graph.Draw("AP") # Draw with axis and points
c.Update()
c.SaveAs("dat_plot.png")

# Keep the canvas open for interactive viewing (optional)
# ROOT.gApplication.Run()