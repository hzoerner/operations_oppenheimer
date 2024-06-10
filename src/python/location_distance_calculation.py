import geopy.distance as distance
import pandas as pd
import os
print(os.getcwd())
from shapely import centroid
from shapely.geometry import Point, Polygon, MultiPolygon
import geopandas as gpd
import shapely
from pyproj import Transformer

admin_map_path = 'operations_oppenheimer/data/vg2500_12-31.utm32s.shape/vg2500/VG2500_KRS.shp'
excel_path = 'operations_oppenheimer/data/ExtendedNuclearData.xlsx'
distance_sheet = "Transport"
final_storage_sheet = "End Storage"
cost_factor = 75
nuclear_factor = 100
cost_factor *= nuclear_factor

# read admin area df
admin_area_df = gpd.read_file(admin_map_path)

# remove cities
admin_area_df = admin_area_df[admin_area_df.BEZ.isin(['Kreis', 'Landkreis'])]
# Define the UTM32 and WGS84 projections
utm32_crs = "EPSG:32632"  # UTM zone 32N
wgs84_crs = "EPSG:4326"   # WGS 84

# Create a transformer object
transformer = Transformer.from_crs(utm32_crs, wgs84_crs)

# Function to transform geometries
def transform_geometry(geom):
    if geom.geom_type == 'Point':
        x, y = transformer.transform(geom.x, geom.y)
        return Point(x, y)
    elif geom.geom_type == 'Polygon':
        exterior = [transformer.transform(x, y) for x, y in geom.exterior.coords]
        interiors = [[transformer.transform(x, y) for x, y in interior.coords] for interior in geom.interiors]
        return Polygon(exterior, interiors)
    elif geom.geom_type == 'MultiPolygon':
        transformed_polygons = []
        for poly in geom.geoms:
            exterior = [transformer.transform(x, y) for x, y in poly.exterior.coords]
            interiors = [[transformer.transform(x, y) for x, y in interior.coords] for interior in poly.interiors]
            transformed_polygons.append(Polygon(exterior, interiors))
        return MultiPolygon(transformed_polygons)
    else:
        raise TypeError(f"Geometry type '{geom.geom_type}' not supported")

# Apply the transformation to the geometry column
admin_area_df['geometry'] = admin_area_df['geometry'].apply(transform_geometry)

# calculate centroids of geometries
admin_area_df["centroids"] = admin_area_df.geometry.apply(lambda x: centroid(x))

### add coordinates of reactors / ISF
# geo coordinates for each interim storage facility
isf_coords = {"Gorleben": Point(53.033344625111255, 11.341597461918763),
              "Ahaus": Point(52.075884172859716, 7.056252098203471),
              "Brokdorf": Point(53.85053568702555, 9.345557443628627),
              "Brunsbüttel": Point(53.891632640377, 9.20050890913132),
              "Biblis": Point(49.70830596704328, 8.411540078513102),
              "Grafenrheinfeld": Point(49.98379398601132, 10.1866276073602),
              "Grohnde": Point(52.03567047221345, 9.41068601909146),
              "Grundemmingen": Point(48.51666979391695, 10.403076045930131),
              "Isar": Point(48.6075745368956, 12.291958334289768),
              "Lingen": Point(52.47033938642687, 7.321564619110349),
              "Krümmel": Point(53.410146129098514, 10.410102309356501),
              "Neckarwestheim": Point(49.0420294948931, 9.173420061294042),
              "Philippsburg": Point(49.25285028023963, 8.441867563151883),
              "Unterweser": Point(53.430121612826966, 8.476324972847728),
              "Lubmin": Point(54.14141668772474, 13.677069821066244)}

# transform regional data in dict format
admin_coords = {k:v for k, v in zip(admin_area_df["GEN"], admin_area_df["centroids"])}

combined_coords = isf_coords | admin_coords

# calculate distances for regular transport
rows_list = []
for f_name, f_point in combined_coords.items():
    f_c = (f_point.x, f_point.y)
    for t_name, t_point in combined_coords.items():
        t_c = (t_point.x, t_point.y)
        rows_list.append({"from": f_name, "to": t_name, "is_possible": 1 if f_name != t_name else 0, "costs": max(375, cost_factor * distance.distance(f_c, t_c).kilometers)})
distance_df = pd.DataFrame(rows_list)

# save distance data to excel file     
with pd.ExcelWriter(excel_path, mode="a", engine="openpyxl", if_sheet_exists="replace", ) as writer:
    distance_df.to_excel(writer, sheet_name=distance_sheet, float_format="%.2f", index=False)


# calculate distances for final storage transport
fsf_coords = {"Saldenberg": Point(48.77298963391,13.35212619416),
              "Lubmin_Endlager": isf_coords["Lubmin"],
              "Messel": Point(49.93831226, 8.749268672)}
storage_coords = isf_coords | admin_coords | fsf_coords
fs_costs = []
for f, f_c in combined_coords.items():
    f_coords = (f_c.x, f_c.y)
    for t, t_c in fsf_coords.items():
        t_coords = (t_c.x, t_c.y)
        fs_costs.append({"from": f, "costs": cost_factor * distance.distance(f_coords, t_coords).kilometers})
fsf_df = pd.DataFrame(fs_costs)

# save distance data to excel file     
with pd.ExcelWriter(excel_path, mode="a", engine="openpyxl", if_sheet_exists="replace", ) as writer:
    fsf_df.to_excel(writer, sheet_name=final_storage_sheet, float_format="%.2f", index=False)


# Update CISF sheet to include all administrative areas
cisf_df = pd.DataFrame()
cisf_df["name"] = list(admin_coords.keys())
cisf_df["capacity"] = 1500
cisf_df["costs"] = 999

with pd.ExcelWriter(excel_path, mode="a", engine="openpyxl", if_sheet_exists="replace", ) as writer:
    cisf_df.to_excel(writer, sheet_name="CISF", float_format="%.2f", index=False)


print("Done.")