import os
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize
from rasterio.windows import Window
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit
from scipy.stats import logistic
from shapely.geometry import Point
import geopandas as gpd

# Répertoires de données
Dir_Obs = "C:/Users/t-perotin/OneDrive - AXA/CAT Modeling/European Windstorms/EU WS ICM - Hazard/2-Reanalysis Selection/"
Dir_Local = "C:/Users/t-perotin/OneDrive - AXA/CAT Modeling/European Windstorms/EU WS ICM - Hazard/6-Statistical Downscaling/"
Dir_Sensitivity = os.path.join(Dir_Local, "SensitivityAnalysis/")
Dir_Predictors = os.path.join(Dir_Local, "Predictors/")

# Chargement des observations
MetOfficeObs = pd.read_csv(os.path.join(Dir_Obs, "FinalMetOfficeObs.csv"))
MetOfficeObs['time'] = pd.to_datetime(MetOfficeObs['time'].str[:13], format="%Y-%m-%dT%H", utc=True)

# Chargement des rasters
rasters = {}
for var in ["Bare", "Crop", "Grass", "Moss", "Shrub", "Snow", "Topo", "Tree", "Urban"]:
    rasters[var] = rasterio.open(os.path.join(Dir_Predictors, f"{var}.01.tif"))

# Calcul des laplaciens et des métriques
def compute_laplacian(data, kernel_size=3):
    kernel = np.array([[0, -1, 0], [-1, 4, -1], [0, -1, 0]])
    return gaussian_filter(data, kernel_size, mode='constant')

# Exemple : Laplacian et terrain roughness
with rasterio.open(os.path.join(Dir_Predictors, "Topo.01.tif")) as src:
    Topo = src.read(1)
    Lapl = compute_laplacian(Topo)
    Roughness = np.gradient(Topo)

# Extraction des prédicteurs locaux
Platforms = MetOfficeObs.groupby('platform_id').agg(
    NObs=('platform_id', 'size'),
    longitude=('longitude', 'mean'),
    latitude=('latitude', 'mean'),
    altitude=('altitude', 'mean'),
    StationStart=('time', 'min'),
    StationEnd=('time', 'max')
).reset_index()

# Création d'un GeoDataFrame
geometry = [Point(xy) for xy in zip(Platforms['longitude'], Platforms['latitude'])]
Platforms = gpd.GeoDataFrame(Platforms, geometry=geometry)
Platforms.set_crs(epsg=4326, inplace=True)

# Extraction des valeurs raster pour les stations
def extract_raster_values(gdf, raster, var_name):
    coords = [(point.x, point.y) for point in gdf.geometry]
    values = [next(raster.sample([coord])) for coord in coords]
    gdf[var_name] = values
    return gdf

for var, raster in rasters.items():
    Platforms = extract_raster_values(Platforms, raster, var)

# Ajustement logistique
def logistic_function(x, a, b, c, d):
    return a / (1 + np.exp(-b * (x - c))) + d

popt, _ = curve_fit(logistic_function, Platforms['Lapl'], Platforms['NObs'])
Platforms['PredictDMU2'] = logistic_function(Platforms['Lapl'], *popt)

# Écriture des rasters ajustés
adjusted_raster_path = os.path.join(Dir_Local, "AdjustedDMU.tif")
with rasterio.open(os.path.join(Dir_Local, "Topo.01.tif")) as src:
    adjusted = logistic_function(src.read(1), *popt)
    with rasterio.open(adjusted_raster_path, 'w', **src.meta) as dst:
        dst.write(adjusted, 1)

# Validation et visualisation
plt.imshow(adjusted, cmap='viridis')
plt.colorbar(label='Adjusted DMU')
plt.title('Adjusted DMU Visualization')
plt.show()

# Enregistrement des résultats
Platforms.to_file(os.path.join(Dir_Local, "Platforms_Enriched.geojson"), driver="GeoJSON")
