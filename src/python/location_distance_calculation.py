import geopy.distance as distance
import pandas as pd

excel_path = 'operations_oppenheimer/ExtendedNuclearData.xlsx'
distance_sheet = "Transport"

# geo coordinates for each interim storage facility
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
              "Unterweser": (53.430121612826966, 8.476324972847728)}

# geo coordinates for hot cells
hot_cell_coords = {"Hot Cell 1": isf_coords["Gorleben"],
                   "Hot Cell 2": isf_coords["Ahaus"]}

#TODO: define reasonale CISF locations

cisf_names = ['Oberdachstetten', 'Kaiseresch', 'Stade', 'Hann. Münden', 'Crivitz', 'Gardelegen', 'Arnstadt', 'Kuchen']
cisf_coords = [(49.41887922, 10.42371430), (50.22679584,7.13714757), (53.59721036,9.45195604), 
(51.41175438, 9.65999936), (53.57227292, 11.64393809), (52.51810692, 11.40461218), (50.83101325, 10.95047937), (48.64054772, 9.802936604)]
cisf_coords = {k:v for k, v in zip(cisf_names, cisf_coords)}


fsf_coords = {"Endlager": (...,...)}

combined_coords = isf_coords | hot_cell_coords | cisf_coords

# calculate distances 
rows_list = []
for f, f_c in combined_coords.items():
    for t, t_c in combined_coords.items():
        rows_list.append({"from": f, "to": t, "is_possible": 1 if f != t else 0, "costs": distance.distance(f_c, t_c).kilometers})
distance_df = pd.DataFrame(rows_list)

# save distance data to excel file     
with pd.ExcelWriter(excel_path, mode="a", engine="openpyxl", if_sheet_exists="replace", ) as writer:
    distance_df.to_excel(writer, sheet_name=distance_sheet, float_format="%.2f", index=False)
#distance_df.head()
print("Done.")