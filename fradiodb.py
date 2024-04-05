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
    remove (dbfilename)
    with sqlite3.connect(dbfilename) as conn:
        cur = conn.cursor()
        for myzipfile in glob(dirpath + sep + "*etalab*.zip"): # [dirpath + sep + x for x in listdir(dirpath) if x.endswith('.zip')]:
            with zipfile.ZipFile(myzipfile) as zFile:
                for csvfile in zFile.infolist():
                    print("importing " + csvfile.filename)
                    tablename = splitext(basename(csvfile.filename))[0].lower()
                    tablename = tablename.replace('sup_', 'anfr_id_' if tablename in ('sup_exploitant','sup_nature','sup_proprietaire','sup_type_antenne') else 'anfr_')
                    pk = {'anfr_support': 'SUP_ID',
                          'anfr_antenne': 'AER_ID',
                          'anfr_station': 'DEM_NM_COMSIS',
                          'anfr_emetteur': 'EMR_ID',
                          'anfr_bande': 'BANEMR_ID',
                          'anfr_id_exploitant': 'ADM_ID',
                          'anfr_id_nature' : 'NAT_ID',
                          'anfr_id_proprietaire' : 'TPO_ID',
                          'anfr_id_type_antenne' : 'TAE_ID'
                        }
                    #if tablename!='anfr_emetteur': continue
                    if(tablename == 'anfr_support'):
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", encoding='iso8859-1', dtype={"STA_NM_ANFR": str, "NAT_ID": 'Int64', "TPO_ID": 'Int64', 'ADR_NM_CP': str}) # "COM_CD_INSEE": 'Int64'} #index_col=pk[tablename],
                        df['COR_NB_LAT'] = round(((df.COR_CD_NS_LAT=='N')*2-1) * (df.COR_NB_DG_LAT + df.COR_NB_MN_LAT/60 + df.COR_NB_SC_LAT/3600), 4)
                        df['COR_NB_LON'] = round(((df.COR_CD_EW_LON=='E')*2-1) * (df.COR_NB_DG_LON + df.COR_NB_MN_LON/60 + df.COR_NB_SC_LON/3600), 4)
                        df['COR_DMS'] = (df.COR_NB_DG_LAT.map(str) + '°' + df.COR_NB_MN_LAT.map(str) + "'" + df.COR_NB_SC_LAT.map(str) + '"' + df.COR_CD_NS_LAT + ' ' +
                                         df.COR_NB_DG_LON.map(str) + "°" + df.COR_NB_MN_LON.map(str) + "'" + df.COR_NB_SC_LON.map(str) + '"' + df.COR_CD_EW_LON)
                        df['ADR_LB_ADD'] = df.ADR_LB_LIEU + df.ADR_LB_ADD1.apply(lambda k: ', ' + k if pd.notna(k) else '') + \
                                                            df.ADR_LB_ADD2.apply(lambda k: ', ' + k if pd.notna(k) else '') + \
                                                            df.ADR_LB_ADD3.apply(lambda k: ', ' + k if pd.notna(k) else '')
                        df[df.ADR_LB_ADD == ''] = None
                        del df['COR_CD_NS_LAT'], df['COR_NB_DG_LAT'], df['COR_NB_MN_LAT'], df['COR_NB_SC_LAT'], df['COR_CD_EW_LON'], df['COR_NB_DG_LON'], df['COR_NB_MN_LON'], df['COR_NB_SC_LON'], df['ADR_LB_ADD1'], df['ADR_LB_ADD2'], df['ADR_LB_ADD3'], df['ADR_LB_LIEU']

                        # SUP_ID is not unique in original ANFR data: one support may host several stations, but also (for historical reasons) one station may be declared with several supports
                        df_stasup = df[['STA_NM_ANFR','SUP_ID']]
                        #df_stasup['STASUP_ID'] = df_stasup.index
                        #df_stasup.set_index('STASUP_ID', inplace=True)
                        df_stasup.to_sql('anfr_stasup', conn, if_exists='replace', index=False)

                        dfsites = df.groupby('COR_DMS', as_index=False)[['COR_DMS','COR_NB_LAT','COR_NB_LON','ADR_LB_ADD','ADR_NM_CP','COM_CD_INSEE']].first()
                        dfsites.to_sql('anfr_site', conn, if_exists='replace', index_label='SITE_ID', dtype={'SITE_ID': 'INTEGER primary key'}) #index=True,
                        dfsites['SITE_ID'] = dfsites.index # FIXME: +1 ?
                        dfsites.set_index('COR_DMS', inplace=True)
                        df['SITE_ID'] = dfsites.loc[df.COR_DMS].SITE_ID.to_numpy()
                        #dfsites['COR_DMS'] = dfsites.index
                        #dfsites.set_index('SITE_ID', inplace=True)

                        dfgroup = df.groupby('SUP_ID', as_index=False)
                        df = dfgroup[['NAT_ID', 'SUP_NM_HAUT', 'TPO_ID', 'SITE_ID']].first()
                        #df['STA_NM_ANFR_list'] = dfgroup['STA_NM_ANFR'].agg(','.join)['STA_NM_ANFR']
                        df.set_index('SUP_ID', inplace=True)
                        #print(dfgroup[dfgroup[['NAT_ID', 'SUP_NM_HAUT', 'TPO_ID', 'COR_DMS']].nunique().ne(1)])
                        #df = dfgroup.agg(','.join)
                    elif(tablename == 'anfr_station'):
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", dtype={"STA_NM_ANFR": str}, parse_dates=['DTE_IMPLANTATION', 'DTE_MODIF'], index_col=pk[tablename], dayfirst=True) #date_format='DD/MM/YYYY', infer_datetime_format=True
                        df['DTE_EN_SERVICE'] = pd.to_datetime(df.DTE_EN_SERVICE, dayfirst=True, errors='coerce') # FIXME: do it more elegantly
                    elif(tablename == 'anfr_emetteur'):
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", dtype={"STA_NM_ANFR": str}, parse_dates=['EMR_DT_SERVICE'], index_col=pk[tablename], dayfirst=True) #infer_datetime_format=True,
                        # systemes = df.groupby('EMR_LB_SYSTEME').size().keys()
                        # df['SYS_ID']=df.EMR_LB_SYSTEME.apply(lambda k: systemes.get_loc(k)+1 if pd.notna(k) else None).astype(pd.Int64Dtype())
                        # del df['EMR_LB_SYSTEME'], df['STA_NM_ANFR'] # FIXME: need some SUP_ID
                        # dfsys=pd.DataFrame(systemes)
                        # dfsys['SYS_ID'] = dfsys.index + 1
                        # dfsys = dfsys[['SYS_ID', 'EMR_LB_SYSTEME']]
                        dfsys = df[['EMR_LB_SYSTEME']].drop_duplicates().sort_values(by='EMR_LB_SYSTEME').reset_index(drop=True)
                        dfsys.to_sql('anfr_id_systeme', conn, if_exists='replace', index_label='SYS_ID', dtype={'SYS_ID': 'INTEGER primary key'}) #index=True,
                        dfsys['SYS_ID'] = dfsys.index
                        dfsys.set_index('EMR_LB_SYSTEME', inplace=True)
                        df['SYS_ID']=dfsys.loc[df.EMR_LB_SYSTEME].SYS_ID.to_numpy()
                        del df['EMR_LB_SYSTEME'], df['STA_NM_ANFR'] # FIXME: need some SUP_ID
                    elif tablename=='anfr_bande':
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", dtype={"STA_NM_ANFR": str}) #, index_col=pk[tablename]
                        df['unit'] = 0
                        df.loc[df.BAN_FG_UNITE=='K', 'unit'] = 1   #1e3
                        df.loc[df.BAN_FG_UNITE=='M', 'unit'] = 1e3 #1e6
                        df.loc[df.BAN_FG_UNITE=='G', 'unit'] = 1e6 #1e9
                        df['BAN_NB_F_DEB'] = round(df.BAN_NB_F_DEB * df.unit).astype(pd.Int64Dtype())
                        df['BAN_NB_F_FIN'] = round(df.BAN_NB_F_FIN * df.unit).astype(pd.Int64Dtype())
                        df.rename(columns={'BAN_ID': 'BANEMR_ID'}, inplace=True)
                        df.set_index('BANEMR_ID', inplace=True)
                        del df['STA_NM_ANFR'], df['BAN_FG_UNITE'], df['unit']
                        dfband = df[['BAN_NB_F_DEB', 'BAN_NB_F_FIN']].drop_duplicates().sort_values(by='BAN_NB_F_DEB').reset_index(drop=True)
                        dfband.to_sql('anfr_id_band', conn, if_exists='replace', index_label='BAN_ID', dtype={'BAN_ID': 'INTEGER primary key'})
                        dfband['BANSTR'] = dfband['BAN_NB_F_DEB'].astype(str) + '_' + dfband['BAN_NB_F_FIN'].astype(str)
                        dfband['BAN_ID'] = dfband.index
                        dfband.set_index('BANSTR', inplace=True)
                        df['BANSTR'] = df['BAN_NB_F_DEB'].astype(str) + '_' + df['BAN_NB_F_FIN'].astype(str)
                        df['BAN_ID'] = dfband.loc[df.BANSTR].BAN_ID.to_numpy()
                        del df['BANSTR'], df['BAN_NB_F_DEB'], df['BAN_NB_F_FIN']
                    elif tablename=='anfr_antenne':
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", dtype={"STA_NM_ANFR": str}) #, index_col=pk[tablename]
                        # AER_ID is not unique in original ANFR data
                        df_staant = df[['STA_NM_ANFR','AER_ID']]
                        df_staant.to_sql('anfr_staant', conn, if_exists='replace', index=False)
                        dfgroup = df.groupby('AER_ID', as_index=False)
                        df = dfgroup[['TAE_ID', 'AER_NB_DIMENSION', 'AER_FG_RAYON', 'AER_NB_AZIMUT', 'AER_NB_ALT_BAS', 'SUP_ID']].first()
                        #df['STA_NM_ANFR_list'] = dfgroup['STA_NM_ANFR'].agg(','.join)['STA_NM_ANFR']
                        df.set_index('AER_ID', inplace=True)
                    else: # id_*
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", dtype={"STA_NM_ANFR": str}, index_col=pk[tablename])

                    cur.execute("drop table if exists " + tablename)
                    df.columns = list(map(lambda x: x.lower(), df.columns))
                    df.to_sql(tablename, conn, if_exists='replace', index_label=pk[tablename], dtype={pk[tablename]: 'INTEGER primary key'})
    create_views(dbfilename)

def create_views(dbfilename):
    with sqlite3.connect(dbfilename) as conn:
        cur = conn.cursor()
        #cur.execute("create view v_support select anfr_support.*, group_concat(sta_nm_anfr,',') sta_nm_anfr_list from anfr_support left outer join anfr_stasup on anfr_stasup.sup_id=anfr_support.SUP_ID group by anfr_support.sup_id")
        cur.executescript("""
                        drop index if exists ix_anfr_stasup_STA_NM_ANFR; create index ix_anfr_stasup_STA_NM_ANFR on anfr_stasup(STA_NM_ANFR);
                        drop index if exists ix_anfr_stasup_SUP_ID; create index ix_anfr_stasup_SUP_ID on anfr_stasup(SUP_ID);
                        drop index if exists ix_anfr_staant_STA_NM_ANFR; create index ix_anfr_staant_STA_NM_ANFR on anfr_staant(STA_NM_ANFR);
                        drop index if exists ix_anfr_staant_AER_ID; create index ix_anfr_staant_AER_ID on anfr_staant(AER_ID);

                        create view v_support as
                        select anfr_support.sup_id, anfr_support.site_id, cor_dms, sup_nm_haut, nat_lb_nom, tpo_lb, sta_nm_anfr_count,sta_nm_anfr_list, aer_id_count,aer_id_list,tech_count,tech_list
                        from anfr_support
                        inner join (
                            select anfr_support.sup_id as sid,count(sta_nm_anfr) sta_nm_anfr_count, group_concat(sta_nm_anfr) sta_nm_anfr_list
                            from anfr_support inner join anfr_stasup on anfr_stasup.sup_id=anfr_support.sup_id
                            group by anfr_support.sup_id
                        ) on sid = anfr_support.sup_id
                        inner join (
                            select anfr_support.sup_id as sid2, count(distinct anfr_antenne.aer_id) aer_id_count, group_concat(distinct anfr_antenne.aer_id) aer_id_list,count(distinct emr_lb_systeme) tech_count, group_concat(distinct emr_lb_systeme) tech_list
                            from anfr_support inner join anfr_antenne on anfr_antenne.sup_id=anfr_support.sup_id
                            inner join anfr_emetteur on anfr_emetteur.aer_id=anfr_antenne.aer_id inner join anfr_id_systeme on anfr_emetteur.sys_id=anfr_id_systeme.sys_id
                            group by anfr_support.sup_id
                        ) on sid2 = anfr_support.sup_id
                        inner join anfr_site on anfr_site.site_id=anfr_support.site_id
                        inner join anfr_id_nature on anfr_id_nature.nat_id=anfr_support.nat_id
                        inner join anfr_id_proprietaire on anfr_id_proprietaire.tpo_id=anfr_support.tpo_id
                        group by anfr_support.sup_id;

                        create view v_antenne as
                        select anfr_antenne.aer_id, sup_id, aer_nb_dimension, aer_fg_rayon, aer_nb_azimut, aer_nb_alt_bas,sta_nm_anfr_count, sta_nm_anfr_list, tae_lb, emr_id_count,emr_id_list,tech_count,tech_list
                        from anfr_antenne
                        inner join (
                            select anfr_antenne.aer_id as aid,count(sta_nm_anfr) sta_nm_anfr_count, group_concat(sta_nm_anfr) sta_nm_anfr_list
                            from anfr_antenne inner join anfr_staant on anfr_staant.aer_id=anfr_antenne.aer_id
                            group by anfr_antenne.aer_id
                        ) on aid = anfr_antenne.aer_id
                        inner join (
                            select anfr_antenne.aer_id as aid2,count(emr_id) emr_id_count, group_concat(emr_id) emr_id_list,
                            count(distinct emr_lb_systeme) tech_count, group_concat(distinct emr_lb_systeme) tech_list
                            from anfr_antenne inner join anfr_emetteur on anfr_emetteur.aer_id=anfr_antenne.aer_id
                            inner join anfr_id_systeme on anfr_emetteur.sys_id=anfr_id_systeme.sys_id
                            group by anfr_antenne.aer_id
                        ) on aid2 = anfr_antenne.aer_id
                        inner join anfr_id_type_antenne on anfr_id_type_antenne.tae_id=anfr_antenne.tae_id
                        group by anfr_antenne.aer_id; """)

def create_imtsectors_table(dbfilename, systype="IMT", bbox="lat<52 and lat>42 and lon<9 and lon>-5", dotable=False):
    tabletype = "table" if dotable else "view"
    with sqlite3.connect(dbfilename) as conn:
        cur = conn.cursor()
        condlist = []
        if bbox:
            condlist.append(bbox)
        if systype:
            if systype=="IMT":
                syslist = [k[0] for k in cur.execute("select emr_lb_systeme from anfr_id_systeme where emr_lb_systeme like 'GSM %' or emr_lb_systeme like 'UMTS %' or emr_lb_systeme like 'LTE %' or emr_lb_systeme like '5G %' ").fetchall()]
                sysnames = "('" + "', '".join(syslist) + "', 'BLR 3 GHZ', 'BLR 3 GHz', 'BLR LTE 3500')"
            elif systype=="BLR": sysnames = "('BLR 3 GHZ', 'BLR 3 GHz', 'BLR LTE 3500')"
            elif systype=="FH": sysnames = "('FH', 'FH ABI')"
            elif systype.startswith('('): sysnames = systype
            condlist.append("emr_lb_systeme in " + sysnames)
        sqlwhere = "where " + " and ".join(condlist)

        print(f'Generating {tabletype} "sectors"')
        tablename="gen_sectors"
        if dotable:
            cur.executescript('''
            drop index if exists anfr_lat_idx;
            drop index if exists anfr_lon_idx;
            drop index if exists anfr_comsis_idx;
            drop index if exists anfr_insee_idx;
            drop index if exists anfr_freq_idx;
            drop index if exists anfr_op_idx;
            ''')
        cur.executescript(f'''
            drop {tabletype} if exists {tablename};
            create {tabletype} {tablename} as
            select anfr_emetteur.emr_id sectorid, dem_nm_comsis comsis,
            anfr_support.sup_id siteid, anfr_station.sta_nm_anfr sta_anfr,
            aer_nb_alt_bas h, sup_nm_haut h2, (ban_nb_f_deb+ban_nb_f_fin)/2e6 freq_mhz,
            (ban_nb_f_fin-ban_nb_f_deb)/1e6 bw_mhz,
            aer_nb_azimut azh, com_cd_insee insee, cor_dms as dms,
            cor_nb_lat as lat, cor_nb_lon as lon, emr_lb_systeme tech,
            adm_lb_nom op, tpo_lb op2, nat_lb_nom as sup_type
            from anfr_bande
            inner join anfr_emetteur on anfr_bande.emr_id=anfr_emetteur.emr_id
            inner join anfr_antenne on anfr_emetteur.aer_id=anfr_antenne.aer_id
            inner join anfr_support on anfr_antenne.sup_id=anfr_support.sup_id
            inner join anfr_station on anfr_antenne.sta_nm_anfr=anfr_station.sta_nm_anfr
            inner join anfr_id_exploitant on anfr_station.adm_id=anfr_id_exploitant.adm_id
            inner join anfr_id_proprietaire on anfr_support.tpo_id=anfr_id_proprietaire.tpo_id
            inner join anfr_id_nature on anfr_support.nat_id=anfr_id_nature.nat_id
            -- inner join anfr_id_systeme on anfr_emetteur.sys_id=anfr_id_systeme.sys_id
            {sqlwhere}
            group by sectorid;
        ''') # null as id, null as pow, null as tilt, null as aer, 0 alt, ban_nb_f_deb fmin,ban_nb_f_fin fmax, # inner join cities on com_cd_insee=cities.insee
        if dotable:
            cur.executescript(f'''
            create index anfr_lat_idx on {tablename}(lat);
            create index anfr_lon_idx on {tablename}(lon);
            create index anfr_comsis_idx on {tablename}(comsis);
            create index anfr_insee_idx on {tablename}(insee);
            create index anfr_freq_idx on {tablename}(freq_mhz);
            create index anfr_op_idx on {tablename}(op);
            ''')

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
