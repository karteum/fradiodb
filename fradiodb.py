#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
@author: Adrien DEMAREZ
"""

import sqlite3
from glob import glob
from os import sep,makedirs,remove
from os.path import splitext,basename,exists
import zipfile
import pandas as pd
from numpy import round
import argparse
import urllib.request as req
import sys

# Schema diagram: d2 ./schema.d2 -s -t 100
# 'NC': ('New Caledonia', (164.029605748, -22.3999760881, 167.120011428, -20.1056458473)),
# lon=165E, lat=21S => X=-1358673, Y=-534964

#def sanitycheck(dbfilename):
    # select * from anfr_support where sup_id not in (select sup_id from anfr_antenne)
    # select * from anfr_emetteur where aer_id not in (select aer_id from anfr_antenne)

def download_data(dirpath='anfr'):
    # From https://www.data.gouv.fr/fr/datasets/donnees-sur-les-installations-radioelectriques-de-plus-de-5-watts-1/
    if not exists(dirpath):
        makedirs(dirpath)
    URLS = {
        "anfr_stations.zip": "https://www.data.gouv.fr/fr/datasets/r/71ba9313-6610-47d7-a5b7-ffaf2fc2427b",
        "anfr_stations_ids.zip": "https://www.data.gouv.fr/fr/datasets/r/dbf19e30-f750-4b25-9dd4-ace5e7d266bb"
        #"densites.xlsx": "https://www.insee.fr/fr/statistiques/fichier/6439600/grille_densite_7_niveaux_2023.xlsx" # https://www.insee.fr/fr/information/6439600
    }
    for k,v in URLS.items():
        print(f"Downloading {k}")
        req.urlretrieve(v, filename=dirpath+sep+k)
    print('Download OK')

def coalesce_freqs(df, idfieldname="transmitter_id"):
    # Coalesce overlapping entries of frequency ranges
    emrbands = {}
    kk=0
    df.sort_values(by='fmin_kHz', inplace=True) # FIXME: needed ?
    for row in df.iterrows():
        #sys.stderr.write(f"\rCoalescing freq entries: {100 * kk // len(df)} %")
        kk+=1
        emr_id,fmin,fmax = row[1]
        if pd.isna(fmin) or pd.isna(fmax): continue
        # small trick: (fmin,fmax) is encoded as (real,imag) to avoid nested lists and associated performance issues
        if not emr_id in emrbands: emrbands[emr_id] = [fmin + 1j*fmax]
        else:
            gbands = emrbands[emr_id]
            found=False
            for k in range(len(gbands)):
                gmin = int(gbands[k].real)
                gmax = int(gbands[k].imag)
                if fmin>=gmin and fmin<=gmax and fmax>=gmin and fmax<=gmax: found=True ; break # Already covered
                elif fmin>=gmin and fmin<=gmax: gbands[k] = gmin+1j*fmax ; found=True ; break # Extend upwards
                elif fmax>=gmin and fmax<=gmax: gbands[k] = fmin+1j*gmax ; found=True ; break # Extend downwards
            if found==False:
                emrbands[emr_id].append(fmin + 1j*fmax)
    res = []
    for k in emrbands.keys():
        for fc in emrbands[k]:
            res.append([k , int(fc.real), int(fc.imag)])
    return pd.DataFrame(data=res, columns=[idfieldname,"fmin_kHz","fmax_kHz"])

def fix_insee_postcode(df):
    import geopandas as gpd
    CRS_WGS84=4326
    gdf_data = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat), crs=CRS_WGS84)
    shapefile_metro = gpd.read_file("shp/communes-20220101.shp")
    gdf_metro = gpd.sjoin(gdf_data, shapefile_metro)
    #gdf_communes_com = gpd.read_file("shp/communes-com-20180101.shp")
    shapefile_nc = gpd.read_file("shp/communes-nc.fgb")
    gdf_nc = gpd.sjoin(gdf_data, shapefile_nc)
    #df['city'] = pd.concat([gdf_metro.nom, gdf_nc.nom])
    #df['insee_correct'] = pd.concat([gdf_metro.insee, gdf_nc.inseecode])
    df['insee_correct'] = gdf_metro['inseecode']

    #shapefile_metro2 = gpd.read_file("shp/COMMUNE_FRMETDROM.shp")
    #gdf_metro2 = gpd.sjoin(gdf_data, shapefile_metro2)
    #df['insee_correct2'] = gdf_metro2['INSEE_COM']
    #return gdf_data, shapefile_metro, gdf_metro, gdf_nc
    return df

def dfdiff(df1, df2):
    return pd.concat([df1,df2]).drop_duplicates(keep=False)

def import_anfr_zip(dbfilename, dirpath='anfr', coalesce=False):
    """Import data from zipped files from data.gouv.fr into a local SQLite DB, with some refinements (e.g. convert DMS coordinates to linear)"""
    if exists(dbfilename):
        remove (dbfilename)
    with sqlite3.connect(dbfilename) as conn:
        cur = conn.cursor()
        df_bandgroup_emr = None
        df_transmitters = None
        for myzipfile in glob(dirpath + sep + "*anfr*.zip"): # [dirpath + sep + x for x in listdir(dirpath) if x.endswith('.zip')]:
            with zipfile.ZipFile(myzipfile) as zFile:
                for csvfile in zFile.infolist():
                    table_rename = {'sup_exploitant': 'id_operators',
                                    'sup_nature': 'id_support_types',
                                    'sup_proprietaire': 'id_support_owners',
                                    'sup_type_antenne': 'id_antenna_types',
                                    'sup_bande': 'bands',
                                    'sup_antenne': 'antennas',
                                    'sup_support': 'supports',
                                    'sup_station': 'stations',
                                    'sup_emetteur': 'transmitters' }
                    tablename = table_rename[splitext(basename(csvfile.filename))[0].lower()]
                    #if not tablename in ('bands', 'transmitters'): continue
                    print("importing " + csvfile.filename)
                    if tablename == 'supports':
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", encoding='iso8859-1', dtype={"STA_NM_ANFR": str, "NAT_ID": 'Int64', "TPO_ID": 'Int64', 'COM_CD_INSEE': str, 'ADR_NM_CP': 'Int64'}) # "COM_CD_INSEE": 'Int64'} #index_col=pk[tablename],index_col="SUP_ID",
                        df.rename(columns={"ADR_NM_CP": "postcode", "COM_CD_INSEE": "inseecode", "NAT_ID": "suptype_id", "SUP_NM_HAUT": "sup_height", "TPO_ID": "owner_id", 'STA_NM_ANFR': "station_name", "SUP_ID": "sup_id"}, inplace=True)
                        df['lat'] = round(((df.COR_CD_NS_LAT=='N')*2-1) * (df.COR_NB_DG_LAT + df.COR_NB_MN_LAT/60 + df.COR_NB_SC_LAT/3600), 4)
                        df['lon'] = round(((df.COR_CD_EW_LON=='E')*2-1) * (df.COR_NB_DG_LON + df.COR_NB_MN_LON/60 + df.COR_NB_SC_LON/3600), 4)
                        df['dms'] = (df.COR_NB_DG_LAT.astype(str) + '°' + df.COR_NB_MN_LAT.astype(str) + "'" + df.COR_NB_SC_LAT.astype(str) + '"' + df.COR_CD_NS_LAT + ' ' +
                                     df.COR_NB_DG_LON.astype(str) + "°" + df.COR_NB_MN_LON.astype(str) + "'" + df.COR_NB_SC_LON.astype(str) + '"' + df.COR_CD_EW_LON)
                        df['address'] = df.ADR_LB_LIEU.str.cat(df[["ADR_LB_ADD1", "ADR_LB_ADD2", "ADR_LB_ADD3"]], sep=', ', na_rep='¤').str.replace(', ¤', '').str.replace('¤, ', '').str.title()
                        del df['COR_CD_NS_LAT'], df['COR_NB_DG_LAT'], df['COR_NB_MN_LAT'], df['COR_NB_SC_LAT'], df['COR_CD_EW_LON'], df['COR_NB_DG_LON'], df['COR_NB_MN_LON'], df['COR_NB_SC_LON'], df['ADR_LB_ADD1'], df['ADR_LB_ADD2'], df['ADR_LB_ADD3'], df['ADR_LB_LIEU']

                        # SUP_ID is not unique in original ANFR data: one support may host several stations, but also (for historical reasons) one station may be declared with several supports
                        # FIXME: not needed ?
                        #df_stasup = df[['station_name']].copy()
                        #df_stasup.to_sql('legacy_stasup', conn, if_exists='replace', index=False)  # FIXME: stasup useless since info is in anfr_emetteur ?

                        # FIXME: original data have an issue because there are 1338 sites that have near-duplicates (i.e. same coordinates but different "address" informations that would need to be merged)
                        # select group_concat(id),dms,group_concat(address, '__¤__'),group_concat(postcode, '__¤__'),group_concat(inseecode, '__¤__'), count(dms) c from supports_tmp group by dms having c>1
                        #dfsites1 = df[['dms','lat','lon','address','postcode','inseecode']].drop_duplicates()
                        dfsites = df.groupby('dms', as_index=False)[['dms','lat','lon','address','postcode','inseecode']].first()
                        #diff = pd.concat([dfsites,dfsites2]).drop_duplicates(keep=False)
                        #dfsites = df.groupby('dms', as_index=False).agg({'address': ' ¤¤ '.join, 'lat': "min", 'lon': "min", 'postcode': "min", 'inseecode': "min"})
                        dfsites.index.name = "id"
                        dfsites.to_sql('sites', conn, if_exists='replace', index_label='id', dtype={'id': 'INTEGER primary key'}) #index=True,
                        dfsites['site_id'] = dfsites.index # FIXME: +1 ?
                        dfsites.set_index('dms', inplace=True)
                        df['site_id'] = dfsites.loc[df.dms].site_id.values

                        #df1 = df.groupby('id', as_index=False)[['suptype_id', 'sup_height', 'owner_id', 'site_id']].first()
                        #df['STA_NM_ANFR_list'] = dfgroup['STA_NM_ANFR'].agg(','.join)['STA_NM_ANFR']
                        df.set_index('sup_id', inplace=True)
                        df = df[['suptype_id', 'sup_height', 'owner_id', 'site_id']].drop_duplicates()
                    elif tablename == 'stations':
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", dtype={"STA_NM_ANFR": str}, parse_dates=['DTE_IMPLANTATION', 'DTE_MODIF'], index_col="DEM_NM_COMSIS", dayfirst=True) #date_format='DD/MM/YYYY', infer_datetime_format=True
                        df['DTE_EN_SERVICE'] = pd.to_datetime(df.DTE_EN_SERVICE, dayfirst=True, errors='coerce') # FIXME: do it more elegantly
                        df.rename(columns={"STA_NM_ANFR": "station_name", "ADM_ID": "operator_id", "DTE_IMPLANTATION": "date_created", "DTE_MODIF": "date_modified", "DTE_EN_SERVICE": "date_operational"}, inplace=True)
                    elif tablename == 'transmitters':
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", dtype={"STA_NM_ANFR": str}, parse_dates=['EMR_DT_SERVICE'], index_col="EMR_ID", dayfirst=True) #infer_datetime_format=True,
                        df.rename(columns={"STA_NM_ANFR": "station_name", "AER_ID": "antenna_id", "EMR_DT_SERVICE": "date_switchedon", "EMR_LB_SYSTEME": "system"}, inplace=True)
                        # FIXME: one transceiver to mutiple antennas <-> one antenna to multiple transmitters ?
                        dfsys = df[['system']].drop_duplicates().sort_values(by='system').reset_index(drop=True)
                        dfsys.index.name = "id"
                        dfsys.to_sql('id_systems', conn, if_exists='replace', index_label='id', dtype={'id': 'INTEGER primary key'}) #index=True,
                        dfsys['system_id'] = dfsys.index
                        dfsys.set_index('system', inplace=True)
                        df['system_id']=dfsys.loc[df.system].system_id.values # FIXME: replace "NULL" entry in systemes
                        del df['system']
                        if df_bandgroup_emr is not None:
                            df['bandgroup_id'] = df_bandgroup_emr.bandgroup_id.astype('Int64')
                        else: df_transmitters = df
                    elif tablename=='bands':
                        print("(Preprocessing may take a few minutes...)")
                        df_banemr = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", usecols=["EMR_ID","BAN_NB_F_DEB","BAN_NB_F_FIN","BAN_FG_UNITE"]) #, index_col=pk[tablename], dtype={"STA_NM_ANFR": str, 'EMR_ID':'Int64'}
                        df_banemr.rename(columns={"STA_NM_ANFR": "station_name", "BAN_ID": "transmitter_band_id", "EMR_ID": "transmitter_id"}, inplace=True)
                        #del df['STA_NM_ANFR']   # anfr_bande already includes a field EMR_ID, and anfr_emetteur already has the correspondance EMR_ID<->STA_NM_ANFR (where an EMR_ID is associated to one and only one )
                        df_banemr['unit'] = 0
                        df_banemr.loc[df_banemr.BAN_FG_UNITE=='K', 'unit'] = 1   #1e3
                        df_banemr.loc[df_banemr.BAN_FG_UNITE=='M', 'unit'] = 1e3 #1e6
                        df_banemr.loc[df_banemr.BAN_FG_UNITE=='G', 'unit'] = 1e6 #1e9
                        df_banemr['fmin_kHz'] = round(df_banemr.BAN_NB_F_DEB * df_banemr.unit).astype('Int64') #pd.Int64Dtype()
                        df_banemr['fmax_kHz'] = round(df_banemr.BAN_NB_F_FIN * df_banemr.unit).astype('Int64')
                        del df_banemr['BAN_FG_UNITE'], df_banemr['unit'], df_banemr['BAN_NB_F_DEB'], df_banemr['BAN_NB_F_FIN']
                        df_banemr.dropna(inplace=True)

                        #if coalesce:
                        #    df_banemr = coalesce_freqs(df_banemr) # A little bit long to compute, but this removes ~24000 useless/duplicate entries
                        df_bands = df_banemr[['fmin_kHz', 'fmax_kHz']].drop_duplicates().sort_values(by='fmin_kHz').reset_index(drop=True)
                        df_bands['BANSTR'] = df_bands['fmin_kHz'].astype(str) + '-' + df_bands['fmax_kHz'].astype(str)
                        df_bands['band_id'] = df_bands.index
                        df_bands.set_index('BANSTR', inplace=True)
                        df_banemr['BANSTR'] = df_banemr['fmin_kHz'].astype(str) + '-' + df_banemr['fmax_kHz'].astype(str)
                        df_banemr['band_id'] = df_bands.loc[df_banemr.BANSTR].band_id.values
                        df_bands.set_index('band_id', inplace=True)
                        df_bands.index.name = "id"

                        df_banemr['band_id_str'] = df_banemr['band_id'].astype(str)
                        df_bandgroup_emr = pd.DataFrame(df_banemr.sort_values('band_id').groupby('transmitter_id')['band_id_str'].unique().agg(','.join))
                        df_bandgroup = df_bandgroup_emr.drop_duplicates().reset_index(drop=True)
                        df_bandgroup['bandgroup_id'] = df_bandgroup.index ; df_bandgroup.set_index('band_id_str', inplace=True)
                        df_bandgroup_emr['bandgroup_id'] = df_bandgroup.loc[df_bandgroup_emr.band_id_str].bandgroup_id.values
                        df_bandgroup['band_id_str'] = df_bandgroup.index ; df_bandgroup.set_index('bandgroup_id', inplace=True)
                        df_bandgroup = pd.DataFrame(df_bandgroup.band_id_str.str.split(',').explode().astype('Int64'))

                        df_bandgroup['fmin_kHz'] = df_bands.loc[df_bandgroup.band_id_str].fmin_kHz.values
                        df_bandgroup['fmax_kHz'] = df_bands.loc[df_bandgroup.band_id_str].fmax_kHz.values
                        #df_bandgroup.rename(columns={'band_id_str': 'band_id'}, inplace=True)
                        #df_bands.to_sql('id_bands', conn, if_exists='replace', index=True, dtype={"id": 'INTEGER primary key'})
                        del df_bandgroup['band_id_str'], df_bandgroup_emr['band_id_str']
                        df_bandgroup.reset_index(inplace=True)
                        df_bandgroup = coalesce_freqs(df_bandgroup, "bandgroup_id")
                        df_bandgroup.to_sql('bandgroups', conn, if_exists='replace', index=False)
                        if df_transmitters is not None:
                            df_transmitters['bandgroup_id'] = df_bandgroup_emr.bandgroup_id.astype('Int64')
                            df = df_transmitters ; tablename = "transmitters"
                        else: df = None
                            #df_bandgroup_emr.to_sql('transmitters_bandgroups', conn, if_exists='replace', index=True)
                        #df_banemr.to_sql('dbg_transmitters_bands', conn, if_exists='replace', index=False)
                    elif tablename=='antennas':
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", dtype={"STA_NM_ANFR": str}) #, index_col=pk[tablename]
                        df.rename(columns={"AER_ID": "antenna_id", "STA_NM_ANFR": "station_name", "TAE_ID": "anttype_id", "AER_NB_DIMENSION": "dimension", "AER_NB_AZIMUT": "azimuth", "AER_FG_RAYON": "dim_type", "AER_NB_ALT_BAS": "ant_height"}, inplace=True)
                        # AER_ID is not unique in original ANFR data. FIXME: is STAANT needed ?
                        #df_staant = df[['station_name','antenna_id']]
                        #df_staant.to_sql('legacy_staant', conn, if_exists='replace', index=False)

                        #dfgroup = df.groupby('AER_ID', as_index=False)
                        #df = dfgroup[['TAE_ID', 'AER_NB_DIMENSION', 'AER_FG_RAYON', 'AER_NB_AZIMUT', 'AER_NB_ALT_BAS', 'SUP_ID']].first()
                        df.set_index('antenna_id', inplace=True)
                        df = df[["anttype_id", "dimension", "dim_type", "azimuth", "ant_height", "SUP_ID"]].drop_duplicates()
                        #df['STA_NM_ANFR_list'] = dfgroup['STA_NM_ANFR'].agg(','.join)['STA_NM_ANFR']
                    elif tablename=='id_operators':
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", dtype={"STA_NM_ANFR": str}, index_col='ADM_ID')
                        df.rename(columns={"ADM_LB_NOM": "operator"}, inplace=True)
                    elif tablename=='id_support_types':
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", dtype={"STA_NM_ANFR": str}, index_col='NAT_ID')
                        df.rename(columns={"NAT_LB_NOM": "support_type"}, inplace=True)
                    elif tablename=='id_support_owners':
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", dtype={"STA_NM_ANFR": str}, index_col='TPO_ID')
                        df.rename(columns={"TPO_LB": "support_owner"}, inplace=True)
                    elif tablename=='id_antenna_types':
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", dtype={"STA_NM_ANFR": str}, index_col='TAE_ID')
                        df.rename(columns={"TAE_LB": "antenna_type"}, inplace=True)
                    else:
                        print(f"Unknown type {tablename}")

                    if df is not None:
                        cur.execute("drop table if exists " + tablename)
                        df.columns = list(map(lambda x: x.lower(), df.columns))
                        df.index.name = "id"
                        #df.to_sql(tablename, conn, if_exists='replace', index_label=pk[tablename], dtype={pk[tablename]: 'INTEGER primary key'})
                        df.to_sql(tablename, conn, if_exists='replace', index_label="id", dtype={"id": 'INTEGER primary key'})
    #create_views(dbfilename)

def gen_sites(dbfilename, dbfilename2=None):
    if dbfilename2==dbfilename: dbfilename2=None
    with sqlite3.connect(dbfilename) as conn:
        cur = conn.cursor()
        global SYSLIST
        SYSLIST = [k[0] for k in cur.execute('select system from id_systems').fetchall()]
        conn.create_function("mask_low", 1, mask_from_list_low64)
        conn.create_function("mask_high", 1, mask_from_list_high64)
        print("gen_sites")
        attachnewdb=f"attach database '{dbfilename2}' as newdb;" if dbfilename2 is not None else ""
        newdb="newdb." if dbfilename2 is not None else ""
        cur.executescript(f"""{attachnewdb}
                             drop table if exists {newdb}gen_sites ;
                             create table {newdb}gen_sites as
                                select sites.id, sites.dms, lon, lat,
                                count(distinct supports.id) support_count, group_concat(distinct supports.id) support_list, max(sup_height) h_max,
                                count(distinct station_name) sta_count, group_concat(distinct station_name) sta_list,
                                count(distinct antennas.id) ant_count, group_concat(distinct antennas.id) ant_list,
                                count(distinct system) tech_count, group_concat(distinct system) tech_list,
                                count(distinct transmitters.id) tx_count, group_concat(distinct transmitters.id) tx_list,
                                count(distinct bandstr) band_count, group_concat(distinct bandstr) band_list,
                                inseecode, 0 as bitmask1, 0 as bitmask2

                                from sites
                                inner join antennas on antennas.sup_id=supports.id
                                inner join supports on sites.id=supports.site_id
                                inner join transmitters on transmitters.antenna_id=antennas.id
                                inner join id_systems on transmitters.system_id=id_systems.id
                                inner join (select bandgroup_id, (fmin_kHz||'-'||fmax_kHz) bandstr from bandgroups) foo on transmitters.bandgroup_id=foo.bandgroup_id
                                group by sites.id;

                             update {newdb}gen_sites set bitmask1=mask_low(tech_list);
                             update {newdb}gen_sites set bitmask2=mask_high(tech_list);

                             drop table if exists {newdb}gen_sectors ;
                             create table {newdb}gen_sectors as
                                select antennas.id, azimuth, ant_height, lon, lat, antennas.sup_id,
                                sup_height, transmitters.station_name, operator,
                                count(distinct system) tech_count, group_concat(distinct system) tech_list,
                                count(distinct transmitters.id) tx_count, group_concat(distinct transmitters.id) tx_list,
                                count(distinct bandstr) band_count, group_concat(distinct bandstr) band_list,
                                0 as bitmask1, 0 as bitmask2

                                from antennas
                                inner join supports on antennas.sup_id=supports.id
                                inner join sites on sites.id=supports.site_id
                                inner join transmitters on transmitters.antenna_id=antennas.id
                                inner join id_systems on transmitters.system_id=id_systems.id
                                inner join stations on stations.station_name=transmitters.station_name
                                inner join (select bandgroup_id, (fmin_kHz||'-'||fmax_kHz) bandstr from bandgroups) foo on transmitters.bandgroup_id=foo.bandgroup_id
                                inner join id_operators on id_operators.id=stations.operator_id
                                group by antennas.id;

                             update {newdb}gen_sectors set bitmask1=mask_low(tech_list);
                             update {newdb}gen_sectors set bitmask2=mask_high(tech_list);
                          """)

def mask_from_list_low64(strlist, masklist=None):
    if pd.isna(strlist): return None
    if masklist is None: masklist = SYSLIST
    mask = 0
    for entry in strlist.split(','):
        mask |= (1<<masklist.index(entry))
    return mask128_low64(mask)
def mask_from_list_high64(strlist, masklist=None):
    if pd.isna(strlist): return None
    if masklist is None: masklist = SYSLIST
    mask = 0
    for entry in strlist.split(','):
        mask |= (1<<masklist.index(entry))
    return mask128_high64(mask)

def list_from_mask(mask, masklist=None):
    if masklist is None: masklist = SYSLIST
    alist = []
    for bit in range(len(masklist)):
        if mask & (1<<bit):
            alist.append(masklist[bit])
    return ','.join(alist)

def mask128_low64(mask128): return mask128 & ((1<<63)-1)  # 63 because sqlite does not deal well with unsigned 64 ints
def mask128_high64(mask128): return mask128 >> 63
def masks64_to_mask128(mask_low, mask_high): return mask_low | (mask_high << 63)

def tablength(dbfilename):
    """List DB tables and size by decreasing number of entries"""
    with sqlite3.connect(dbfilename) as conn:
        cur = conn.cursor()
        tables = [k[0] for k in cur.execute("select name from sqlite_schema where type='table'").fetchall()] # FROM sqlite_master ?
        #tables = list(zip(*cur.execute("select name from sqlite_schema where type='table'").fetchall()))[0]
        res = {}
        for table in tables:
            #elements = [k[0] for k in cur.execute(f"SELECT count(*) FROM {table}").fetchall()][0]
            #tmp[table] = elements
            length = cur.execute(f"select count(*) from {table}").fetchall()[0][0]
            res[length] = table
    #return pd.Series(tmp).sort_values()
    for k in sorted(res.keys(), reverse=True):
        print(f"{res[k]} : {k}")
    return res

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", "-d", help="Download data from data.gouv.fr", action='store_true', default=False)
    parser.add_argument("--datadir", "-p", help="Data location", default="anfr")
    parser.add_argument("--importdb", "-i", help="Import CSV into DB", default=None)
    parser.add_argument("--gensectorsdb", "-g", help="Create additional table/DB with summarized data", default=None)

    args = parser.parse_args()
    if args.download:
        download_data(args.datadir)
    if args.importdb is not None:
        #import_cities(args.dbfile, dirpath=mydir)
        import_anfr_zip(args.importdb, dirpath=args.datadir)
    if args.gensectorsdb is not None:
        origdb = args.importdb if args.importdb is not None else args.gensectorsdb
        gen_sites(origdb, args.gensectorsdb)
