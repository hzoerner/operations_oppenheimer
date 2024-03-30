import pandas as pd
from shapely.geometry import Point
import geopandas as gpd
from geopandas import GeoDataFrame
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


isf_coords = {"Gorleben": (53.033344625111255, 11.341597461918763),
              "Ahaus": (52.075884172859716, 7.056252098203471),
              "Brokdorf": (53.85053568702555, 9.345557443628627),
              "Brunsbüttel": (53.891632640377, 9.20050890913132),
              "Biblis": (49.70830596704328, 8.411540078513102),
              "Grafenrheinfeld": (49.98379398601132, 10.1866276073602),
              "Grohnde": (52.03567047221345, 9.41068601909146),
              "Grundemmingen": (48.51666979391695, 10.403076045930131),
              "Isar": (48.6075745368956, 12.291958334289768),
              "Lingen": (52.47033938642687, 7.321564619110349),
              "Krümmel": (53.410146129098514, 10.410102309356501),
              "Neckarwestheim": (49.0420294948931, 9.173420061294042),
              "Philippsburg": (49.25285028023963, 8.441867563151883),
              "Unterweser": (53.430121612826966, 8.476324972847728),
              "Lubmin": (54.14141668772474, 13.677069821066244)}

cisf_names = ['Oberdachstetten', 'Kaiseresch', 'Stade', 'Hann. Münden', 'Crivitz', 'Gardelegen', 'Arnstadt', 'Kuchen']
cisf_coords = [(49.41887922, 10.42371430), (50.22679584,7.13714757), (53.59721036,9.45195604), 
(51.41175438, 9.65999936), (53.57227292, 11.64393809), (52.51810692, 11.40461218), (50.83101325, 10.95047937), (48.64054772, 9.802936604)]
cisf_coords = {k:v for k, v in zip(cisf_names, cisf_coords)}

fsf_coords = {"Saldenberg": (48.77298963391,13.35212619416),
              "Gorleben_Endlager": (53.033344625111255, 11.341597461918763),
              "Messel": (49.93831226, 8.749268672)}

coords = isf_coords | cisf_coords | fsf_coords

isf_snf = pd.read_excel("operations_oppenheimer/data/ExtendedNuclearData.xlsx", sheet_name="Reactors")
marker_sizes = [isf_snf[isf_snf.name == n]["snf"].values[0] for n in isf_coords.keys()]
row_list = []
for n, c in coords.items():
    d = {}
    d["name"] = n
    d["latitude"] = c[0]
    d["longitude"] = c[1]
    row_list.append(d)

coords_df = pd.DataFrame(row_list)

# plot data
# generate axis handles
teal_star = mlines.Line2D([], [], color='teal', marker='*', linestyle='None',
                          markersize=10, label='FSF')
tomato_x = mlines.Line2D([], [], color='tomato', marker='x', linestyle='None',
                          markersize=10, label='CISF')
blue_circle = mlines.Line2D([], [], color='deepskyblue', marker='o', linestyle='None',
                          markersize=10, label='ISF')

# read geodata
geometry = [Point(xy) for xy in zip(coords_df['longitude'], coords_df['latitude'])]
gdf = GeoDataFrame(coords_df, geometry=geometry)
ger = gpd.read_file("operations_oppenheimer/data/vg2500_geo84/vg2500_krs.shp")


# plot existing ISF
isf_ax = gdf[gdf.name.isin(isf_coords.keys())].plot(aspect=1.4, ax=ger.plot(color="black", alpha=0.25, label="ISF"),
         c="deepskyblue", 
         markersize=marker_sizes, label="ISF")
isf_ax.set_title("Storage locations")
isf_ax.set_axis_off()
plt.legend(handles=[blue_circle], fancybox=True, loc='best', bbox_to_anchor=(0.5, 0., 0.5, 0.5))
plt.savefig("reactors.png")

# plot CISF
isf_ax = gdf[gdf.name.isin(cisf_coords.keys())].plot(aspect=1.4, ax=isf_ax, marker="x",
         c="tomato", 
         markersize=isf_snf.snf.mean() / 2, label="CISF")
plt.legend(handles=[blue_circle, tomato_x], fancybox=True, loc='best', bbox_to_anchor=(0.5, 0., 0.5, 0.5))
plt.savefig("cisf_newonly.png")

isf_ax = gdf[gdf.name.isin(["Gorleben", "Ahaus"])].plot(aspect=1.4, ax=isf_ax, marker="x",
         c="tomato", 
         markersize=isf_snf.snf.mean() * 3, label="CISF")
plt.legend(handles=[blue_circle, tomato_x], fancybox=True, loc='best', bbox_to_anchor=(0.5, 0., 0.5, 0.5))
plt.savefig("cisf.png")

# plot final storage
isf_ax = gdf[gdf.name.isin(fsf_coords.keys())].plot(aspect=1.4, ax=isf_ax, marker="*", alpha=0.75,
         c="teal", 
         markersize=isf_snf.snf.mean() * 2, label="FSF", legend=True)
plt.legend(handles=[blue_circle, tomato_x, teal_star], fancybox=True, loc='best', bbox_to_anchor=(0.5, 0., 0.5, 0.5))
plt.savefig("complete_locations.png")

#plt.show() 