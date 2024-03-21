import geopy.distance as distance
import pandas as pd

excel_path = 'NuclearData.xlsx'
distance_sheet = "Transport"
# read orginal excel file with costs/distances in Transport sheet
xlsx = pd.ExcelFile(excel_path)
distance_df = pd.read_excel(xlsx, distance_sheet)

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
cisf_coords = {}


combined_coords = isf_coords | hot_cell_coords # NOTE: cisf_coords not included

for f, f_c in combined_coords.items():
    for t, t_c in combined_coords.items():
        if f_c != t_c:
            distance_df.loc[(distance_df["from"] == f) & (distance_df["to"] == t), "is_possible"] = 1
            distance_df.loc[(distance_df["from"] == f) & (distance_df["to"] == t), "costs"] = distance.distance(f_c, t_c).kilometers
        else:
            distance_df.loc[(distance_df["from"] == f) & (distance_df["to"] == t), "is_possible"] = 0
            distance_df.loc[(distance_df["from"] == f) & (distance_df["to"] == t), "costs"] = 0

# save updated data to original excel file            
distance_df.to_excel(excel_path, sheet_name=distance_sheet, float_format="%.2f", index=False)
distance_df.head()