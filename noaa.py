import pandas as pd
import geopandas as gpd
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.parse import urljoin
from tqdm.auto import tqdm

# Constantes
BASE_URL = "https://www.ncei.noaa.gov/pub/data/noaa/isd-lite/{year}/"
FIELDS = ('temp', 'dewtemp', 'pressure', 'winddirection', 'windspeed', 'skycoverage', 'precipitation-1h', 'precipitation-6h')
METADATA_URL = 'https://www.ncei.noaa.gov/pub/data/noaa/isd-history.txt'

def fetch_metadata():
    """Télécharge et traite les métadonnées des stations météo."""
    metadata = pd.read_fwf(METADATA_URL, skiprows=19)
    metadata = metadata.dropna(subset=['LAT', 'LON'])
    metadata = metadata[~((metadata.LON == 0) & (metadata.LAT == 0))]
    metadata['BEGIN'] = pd.to_datetime(metadata['BEGIN'].astype(str), errors='coerce')
    metadata['END'] = pd.to_datetime(metadata['END'].astype(str), errors='coerce')
    metadata = metadata.dropna(subset=['BEGIN', 'END'])
    metadata_gdf = gpd.GeoDataFrame(metadata, geometry=gpd.points_from_xy(metadata.LON, metadata.LAT, crs="EPSG:4326"))
    return metadata_gdf

def filter_stations(metadata, countries, start, end):
    """Filtre les stations par pays et période de disponibilité."""
    if isinstance(countries, str):
        countries = [countries]
    filtered = metadata[metadata['CTRY'].isin(countries)]
    filtered = filtered[(filtered['BEGIN'] <= pd.to_datetime(start)) & (filtered['END'] >= pd.to_datetime(end))]
    return filtered['USAF'].unique()

def download_station_data(station_id, years):
    """Télécharge les données d'une station pour les années spécifiées."""
    dfs = []
    for year in years:
        try:
            url = urljoin(BASE_URL.format(year=year), f"{station_id}-99999-{year}.gz")
            df = pd.read_csv(url, sep='\\s+', header=None, na_values=-9999)
            df.columns = ['year', 'month', 'day', 'hour'] + list(FIELDS)
            df[['temp', 'dewtemp', 'pressure', 'windspeed']] /= 10.0
            df['timestamp'] = pd.to_datetime(df[['year', 'month', 'day', 'hour']])
            df = df.set_index('timestamp').drop(columns=['year', 'month', 'day', 'hour'])
            dfs.append(df)
        except Exception:
            continue
    return pd.concat(dfs) if dfs else pd.DataFrame()

def organize_by_field(data_dict):
    """Réorganise les données par champ météo."""
    organized = {field: pd.concat([data[field].rename(station_id) for station_id, data in data_dict.items()], axis=1) 
                 for field in FIELDS}
    return organized

def f(start, end, countries):
    """
    Récupère les données météo ISD-Lite pour une période et des pays donnés, organisées par champ.
    
    Args:
        start (str or datetime): Date de début ('YYYY-MM-DD').
        end (str or datetime): Date de fin ('YYYY-MM-DD').
        countries (str or list of str): Codes pays ISO alpha-2 (e.g., 'US', ['FR', 'DE']).
        
    Returns:
        dict: Données organisées par champ, sous forme de DataFrames.
    """
    # Téléchargement et filtrage des métadonnées
    metadata = fetch_metadata()
    station_ids = filter_stations(metadata, countries, start, end)
    
    # Calcul des années à couvrir
    start, end = pd.to_datetime(start), pd.to_datetime(end)
    years = list(range(start.year, end.year + 1))
    
    # Téléchargement des données pour chaque station
    data_by_station = {}
    with ThreadPoolExecutor(max_workers=6) as executor:
        futures = {executor.submit(download_station_data, station_id, years): station_id for station_id in station_ids}
        for future in tqdm(as_completed(futures), total=len(futures), desc="Downloading data"):
            station_id = futures[future]
            try:
                data = future.result()
                if not data.empty:
                    data_by_station[station_id] = data
            except Exception:
                continue

    # Organisation des données par champ météo
    return organize_by_field(data_by_station)
