#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
@author: Adrien DEMAREZ
"""

import sqlite3
from glob import glob
from os import sep,makedirs
from os.path import splitext,basename,exists
import zipfile
import pandas as pd
from numpy import round
import argparse
import urllib.request as req

def download_data(dirpath='etalab'):
    # From https://www.data.gouv.fr/fr/datasets/donnees-sur-les-installations-radioelectriques-de-plus-de-5-watts-1/
    if not exists(dirpath):
        makedirs(dirpath)
    URLS = {
        "etalab_stations.zip": "https://www.data.gouv.fr/fr/datasets/r/8500e7f6-ffcd-4c09-840e-2eadd414130e",
        "etalab_stations_ids.zip": "https://www.data.gouv.fr/fr/datasets/r/78e112c9-1958-4e58-8762-1568192517e5",
        "densites.xlsx": "https://www.insee.fr/fr/statistiques/fichier/6439600/grille_densite_7_niveaux_2023.xlsx" # https://www.insee.fr/fr/information/6439600
    }
    for k,v in URLS.items():
        print(f"Downloading {k}")
        req.urlretrieve(v, filename=dirpath+sep+k)
    print('Download OK')

def import_cities(dbfilename, dirpath="etalab"):
    """Import density grid on a per-city basis as made by INSEE on https://www.insee.fr/fr/information/6439600"""
    print("Importing city density grid from INSEE")
    xlsfilename=dirpath + sep + "densites.xlsx"
    df=pd.read_excel(xlsfilename, sheet_name='Grille_Densite', skiprows=4, header=0)
    #df.columns = ['insee', 'name', 'dclass', 'region', 'population', 'dfrac1', 'dfrac2', 'dfrac3', 'dfrac4', 'dfrac5', 'dfrac6', 'dfrac7']
    with sqlite3.connect(dbfilename) as conn:
        cur = conn.cursor()
        cur.execute("drop table if exists cities")
        df.to_sql("cities", conn, index=False)

def import_etalab_zip(dbfilename, dirpath='etalab'):
    """Import data from zipped files from data.gouv.fr into a local SQLite DB, with some refinements (e.g. convert DMS coordinates to linear)"""
    with sqlite3.connect(dbfilename) as conn:
        cur = conn.cursor()
        for myzipfile in glob(dirpath + sep + "*etalab*.zip"): # [dirpath + sep + x for x in listdir(dirpath) if x.endswith('.zip')]:
            with zipfile.ZipFile(myzipfile) as zFile:
                for csvfile in zFile.infolist():
                    print("importing " + csvfile.filename)
                    tablename = splitext(basename(csvfile.filename))[0].lower()
                    tablename = tablename.replace('sup_', 'eta_id_' if tablename in ('sup_exploitant','sup_nature','sup_proprietaire','sup_type_antenne') else 'eta_')
                    pk = {'eta_support': 'SUP_ID', # FIXME: not unique
                          'eta_antenne': 'AER_ID', # FIXME: not unique
                          'eta_station': 'DEM_NM_COMSIS',
                          'eta_emetteur': 'EMR_ID',
                          'eta_bande': 'BAN_ID',
                          'eta_id_exploitant': 'ADM_ID',
                          'eta_id_nature' : 'NAT_ID',
                          'eta_id_proprietaire' : 'TPO_ID',
                          'eta_id_type_antenne' : 'TAE_ID'
                        }
                    if(tablename == 'eta_support'):
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", encoding='iso8859-1', index_col=pk[tablename], dtype={"STA_NM_ANFR": object, "NAT_ID": 'Int64', "TPO_ID": 'Int64', 'ADR_NM_CP': object}) # "COM_CD_INSEE": 'Int64'}
                        df['COR_NB_LAT'] = round(((df.COR_CD_NS_LAT=='N')*2-1) * (df.COR_NB_DG_LAT + df.COR_NB_MN_LAT/60 + df.COR_NB_SC_LAT/3600), 4)
                        df['COR_NB_LON'] = round(((df.COR_CD_EW_LON=='E')*2-1) * (df.COR_NB_DG_LON + df.COR_NB_MN_LON/60 + df.COR_NB_SC_LON/3600), 4)
                        df['COR_DMS'] = (df.COR_NB_DG_LAT.map(str) + '°' + df.COR_NB_MN_LAT.map(str) + "'" + df.COR_NB_SC_LAT.map(str) + '"' + df.COR_CD_NS_LAT + ' ' +
                                         df.COR_NB_DG_LON.map(str) + "°" + df.COR_NB_MN_LON.map(str) + "'" + df.COR_NB_SC_LON.map(str) + '"' + df.COR_CD_EW_LON)
                        df['ADR_LB_ADD'] = df.ADR_LB_LIEU + df.ADR_LB_ADD1.apply(lambda k: ' ; ' + k if pd.notna(k) else '') + \
                                                            df.ADR_LB_ADD2.apply(lambda k: ' ; ' + k if pd.notna(k) else '') + \
                                                            df.ADR_LB_ADD3.apply(lambda k: ' ; ' + k if pd.notna(k) else '')
                        df[df.ADR_LB_ADD == ''] = None
                        # Delete useless fields
                        del df['COR_CD_NS_LAT'], df['COR_NB_DG_LAT'], df['COR_NB_MN_LAT'], df['COR_NB_SC_LAT'], df['COR_CD_EW_LON'], df['COR_NB_DG_LON'], df['COR_NB_MN_LON'], df['COR_NB_SC_LON'], df['ADR_LB_ADD1'], df['ADR_LB_ADD2'], df['ADR_LB_ADD3'], df['ADR_LB_LIEU']
                    elif(tablename == 'eta_station'):
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", dtype={"STA_NM_ANFR": object}, parse_dates=['DTE_IMPLANTATION', 'DTE_MODIF', 'DTE_EN_SERVICE'], index_col=pk[tablename], dayfirst=True) #date_format='DD/MM/YYYY', infer_datetime_format=True
                    elif(tablename == 'eta_emetteur'):
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", dtype={"STA_NM_ANFR": object}, parse_dates=['EMR_DT_SERVICE'], index_col=pk[tablename], dayfirst=True) #infer_datetime_format=True,
                        systemes = df.groupby('EMR_LB_SYSTEME').size().keys()
                        df['SYS_ID']=df.EMR_LB_SYSTEME.apply(lambda k: systemes.get_loc(k)+1 if pd.notna(k) else None).astype(pd.Int64Dtype())
                        del df['EMR_LB_SYSTEME'], df['STA_NM_ANFR'] # FIXME: need some SUP_ID
                        dfsys=pd.DataFrame(systemes)
                        dfsys['SYS_ID'] = dfsys.index + 1
                        dfsys = dfsys[['SYS_ID', 'EMR_LB_SYSTEME']]
                        dfsys.to_sql('eta_id_systeme', conn, if_exists='append', index=False, index_label='SYS_ID')
                    elif tablename=='eta_bande':
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", dtype={"STA_NM_ANFR": object}, index_col=pk[tablename])
                        df['unit'] = 0
                        df.loc[df.BAN_FG_UNITE=='K', 'unit'] = 1e3
                        df.loc[df.BAN_FG_UNITE=='M', 'unit'] = 1e6
                        df.loc[df.BAN_FG_UNITE=='G', 'unit'] = 1e9
                        df['BAN_NB_F_DEB'] = round(df.BAN_NB_F_DEB * df.unit).astype(pd.Int64Dtype())
                        df['BAN_NB_F_FIN'] = round(df.BAN_NB_F_FIN * df.unit).astype(pd.Int64Dtype())
                        del df['STA_NM_ANFR'], df['BAN_FG_UNITE'], df['unit']
                    else: # antenne & id_*
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", dtype={"STA_NM_ANFR": object}, index_col=pk[tablename])
                    cur.execute("drop table if exists " + tablename)
                    df.columns = list(map(lambda x: x.lower(), df.columns))
                    if tablename!='eta_antenne' and tablename!='eta_support':
                        df.to_sql(tablename, conn, if_exists='replace', index_label=pk[tablename], dtype={pk[tablename]: 'INTEGER primary key'})
                    else:
                        df.to_sql(tablename, conn, if_exists='replace', index_label=pk[tablename]) # FIXME: solve issues on PK


def create_imtsectors_table(dbfilename, systype="IMT", bbox="lat<52 and lat>42 and lon<9 and lon>-5", dotable=False):
    tabletype = "table" if dotable else "view"
    with sqlite3.connect(dbfilename) as conn:
        cur = conn.cursor()
        condlist = []
        if bbox:
            condlist.append(bbox)
        if systype:
            if systype=="IMT":
                syslist = [k[0] for k in cur.execute("select emr_lb_systeme from eta_id_systeme where emr_lb_systeme like 'GSM %' or emr_lb_systeme like 'UMTS %' or emr_lb_systeme like 'LTE %' or emr_lb_systeme like '5G %' ").fetchall()]
                sysnames = "('" + "', '".join(syslist) + "', 'BLR 3 GHZ', 'BLR 3 GHz', 'BLR LTE 3500')"
            elif systype=="BLR": sysnames = "('BLR 3 GHZ', 'BLR 3 GHz', 'BLR LTE 3500')"
            elif systype=="FH": sysnames = "('FH', 'FH ABI')"
            elif systype.startswith('('): sysnames = systype
            condlist.append("eta_id_systeme.emr_lb_systeme in " + sysnames)
        sqlwhere = "where " + " and ".join(condlist)

        print(f'Generating {tabletype} "sectors"')
        tablename="gen_sectors"
        if dotable:
            cur.executescript('''
            drop index if exists eta_lat_idx;
            drop index if exists eta_lon_idx;
            drop index if exists eta_comsis_idx;
            drop index if exists eta_insee_idx;
            drop index if exists eta_freq_idx;
            drop index if exists eta_op_idx;
            ''')
        cur.executescript(f'''
            drop {tabletype} if exists {tablename};
            create {tabletype} {tablename} as
            select eta_emetteur.emr_id sectorid, dem_nm_comsis comsis,
            eta_support.sup_id siteid, eta_station.sta_nm_anfr sta_anfr,
            aer_nb_alt_bas h, sup_nm_haut h2, (ban_nb_f_deb+ban_nb_f_fin)/2e6 freq_mhz,
            (ban_nb_f_fin-ban_nb_f_deb)/1e6 bw_mhz,
            aer_nb_azimut azh, com_cd_insee insee, cor_dms as dms,
            cor_nb_lat as lat, cor_nb_lon as lon, emr_lb_systeme tech,
            adm_lb_nom op, tpo_lb op2, nat_lb_nom as sup_type
            from eta_bande
            inner join eta_id_exploitant on eta_station.adm_id=eta_id_exploitant.adm_id
            inner join eta_id_proprietaire on eta_support.tpo_id=eta_id_proprietaire.tpo_id
            inner join eta_id_nature on eta_support.nat_id=eta_id_nature.nat_id
            inner join eta_emetteur on eta_bande.emr_id=eta_emetteur.emr_id
            inner join eta_id_systeme on eta_emetteur.sys_id=eta_id_systeme.sys_id
            inner join eta_antenne on eta_antenne.aer_id=eta_emetteur.aer_id
            inner join eta_support on eta_antenne.sta_nm_anfr=eta_support.sta_nm_anfr
            inner join eta_station on eta_antenne.sta_nm_anfr=eta_station.sta_nm_anfr
            {sqlwhere} group by sectorid;
        ''') # null as id, null as pow, null as tilt, null as aer, 0 alt, ban_nb_f_deb fmin,ban_nb_f_fin fmax, # inner join cities on com_cd_insee=cities.insee
        if dotable:
            cur.executescript(f'''
            create index eta_lat_idx on {tablename}(lat);
            create index eta_lon_idx on {tablename}(lon);
            create index eta_comsis_idx on {tablename}(comsis);
            create index eta_insee_idx on {tablename}(insee);
            create index eta_freq_idx on {tablename}(freq_mhz);
            create index eta_op_idx on {tablename}(op);
            ''')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dbfile", help="DB path")
    parser.add_argument("--download", "-d", help="Display files from dirA that are in dirB", action='store_true', default=False)
    parser.add_argument("--usedir", "-u", help="Data location", action='store_true', default=False)
    parser.add_argument("--sectors", "-s", help="Make table 'sectors'", action='store_true', default=False)
    parser.add_argument("--reset", "-r", help="Reset DB", action='store_true', default=False)

    args = parser.parse_args()
    dirpath = args.usedir if args.usedir else "etalab"
    if args.download:
        download_data(dirpath)

    if args.reset:
        #import_cities(args.dbfile, dirpath=mydir)
        import_etalab_zip(args.dbfile, dirpath=dirpath)
    if args.sectors:
        create_imtsectors_table(args.dbfile)
