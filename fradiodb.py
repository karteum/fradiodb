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

#def sanitycheck(dbfilename):
    # select * from anfr_support where sup_id not in (select sup_id from anfr_antenne)
    # select * from anfr_emetteur where aer_id not in (select aer_id from anfr_antenne)

def download_data(dirpath='etalab'):
    # From https://www.data.gouv.fr/fr/datasets/donnees-sur-les-installations-radioelectriques-de-plus-de-5-watts-1/
    if not exists(dirpath):
        makedirs(dirpath)
    URLS = {
        "etalab_stations.zip": "https://www.data.gouv.fr/fr/datasets/r/f03b1594-d0fe-4b2c-8b11-497d014cdcd0",
        "etalab_stations_ids.zip": "https://www.data.gouv.fr/fr/datasets/r/56126eb6-5c22-4de5-a816-e02d5146c7b2",
        "densites.xlsx": "https://www.insee.fr/fr/statistiques/fichier/6439600/grille_densite_7_niveaux_2023.xlsx" # https://www.insee.fr/fr/information/6439600
    }
    for k,v in URLS.items():
        print(f"Downloading {k}")
        req.urlretrieve(v, filename=dirpath+sep+k)
    print('Download OK')

def coalesce_freqs(df):
    # Coalesce overlapping entries of frequency ranges
    emrbands = {}
    kk=0
    for row in df.iterrows():
        sys.stderr.write(f"\rCoalescing freq entries: {100 * kk // len(df)} %")
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
    return pd.DataFrame(data=res, columns=["EMR_ID","BAN_NB_F_DEB","BAN_NB_F_FIN"])

def import_etalab_zip(dbfilename, dirpath='etalab', coalesce=True):
    """Import data from zipped files from data.gouv.fr into a local SQLite DB, with some refinements (e.g. convert DMS coordinates to linear)"""
    if exists(dbfilename):
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
                          'anfr_bande': 'BAN_ID',
                          'anfr_id_exploitant': 'ADM_ID',
                          'anfr_id_nature' : 'NAT_ID',
                          'anfr_id_proprietaire' : 'TPO_ID',
                          'anfr_id_type_antenne' : 'TAE_ID'
                        }
                    #if tablename!='anfr_bande': continue
                    if(tablename == 'anfr_support'):
                        df = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", encoding='iso8859-1', dtype={"STA_NM_ANFR": str, "NAT_ID": 'Int64', "TPO_ID": 'Int64', 'COM_CD_INSEE': str, 'ADR_NM_CP': 'Int64'}) # "COM_CD_INSEE": 'Int64'} #index_col=pk[tablename],
                        df['COR_NB_LAT'] = round(((df.COR_CD_NS_LAT=='N')*2-1) * (df.COR_NB_DG_LAT + df.COR_NB_MN_LAT/60 + df.COR_NB_SC_LAT/3600), 4)
                        df['COR_NB_LON'] = round(((df.COR_CD_EW_LON=='E')*2-1) * (df.COR_NB_DG_LON + df.COR_NB_MN_LON/60 + df.COR_NB_SC_LON/3600), 4)
                        df['COR_DMS'] = (df.COR_NB_DG_LAT.map(str) + '°' + df.COR_NB_MN_LAT.map(str) + "'" + df.COR_NB_SC_LAT.map(str) + '"' + df.COR_CD_NS_LAT + ' ' +
                                         df.COR_NB_DG_LON.map(str) + "°" + df.COR_NB_MN_LON.map(str) + "'" + df.COR_NB_SC_LON.map(str) + '"' + df.COR_CD_EW_LON)
                        df['ADR_LB_ADD'] = df.ADR_LB_LIEU.str.cat(df[["ADR_LB_ADD1", "ADR_LB_ADD2", "ADR_LB_ADD3"]], sep=', ', na_rep='¤').str.replace(', ¤', '').str.replace('¤, ', '')
                        del df['COR_CD_NS_LAT'], df['COR_NB_DG_LAT'], df['COR_NB_MN_LAT'], df['COR_NB_SC_LAT'], df['COR_CD_EW_LON'], df['COR_NB_DG_LON'], df['COR_NB_MN_LON'], df['COR_NB_SC_LON'], df['ADR_LB_ADD1'], df['ADR_LB_ADD2'], df['ADR_LB_ADD3'], df['ADR_LB_LIEU']

                        # SUP_ID is not unique in original ANFR data: one support may host several stations, but also (for historical reasons) one station may be declared with several supports
                        #df_stasup = df[['STA_NM_ANFR','SUP_ID']]
                        #df_stasup.to_sql('gen_stasup', conn, if_exists='replace', index=False)  # FIXME: stasup useless since info is in anfr_emetteur ?

                        # FIXME: check whether TPO_ID ought to be moved in dfsites (are there any sites with multiple supports from different owners ?)
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
                        # FIXME: one transceiver to mutiple antennas <-> one antenna to multiple transmitters ?
                        dfsys = df[['EMR_LB_SYSTEME']].drop_duplicates().sort_values(by='EMR_LB_SYSTEME').reset_index(drop=True)
                        dfsys.to_sql('anfr_id_systeme', conn, if_exists='replace', index_label='SYS_ID', dtype={'SYS_ID': 'INTEGER primary key'}) #index=True,
                        dfsys['SYS_ID'] = dfsys.index
                        dfsys.set_index('EMR_LB_SYSTEME', inplace=True)
                        df['SYS_ID']=dfsys.loc[df.EMR_LB_SYSTEME].SYS_ID.to_numpy() # FIXME: replace "NULL" entry in systemes
                        del df['EMR_LB_SYSTEME']
                    elif tablename=='anfr_bande':
                        df_banemr = pd.read_csv(zFile.open(csvfile), sep=';', decimal = ",", usecols=["EMR_ID","BAN_NB_F_DEB","BAN_NB_F_FIN","BAN_FG_UNITE"]) #, index_col=pk[tablename], dtype={"STA_NM_ANFR": str, 'EMR_ID':'Int64'}
                        #del df['STA_NM_ANFR']   # anfr_bande already includes a field EMR_ID, and anfr_emetteur already has the correspondance EMR_ID<->STA_NM_ANFR (where an EMR_ID is associated to one and only one )
                        df_banemr['unit'] = 0
                        df_banemr.loc[df_banemr.BAN_FG_UNITE=='K', 'unit'] = 1   #1e3
                        df_banemr.loc[df_banemr.BAN_FG_UNITE=='M', 'unit'] = 1e3 #1e6
                        df_banemr.loc[df_banemr.BAN_FG_UNITE=='G', 'unit'] = 1e6 #1e9
                        df_banemr['BAN_NB_F_DEB'] = round(df_banemr.BAN_NB_F_DEB * df_banemr.unit).astype(pd.Int64Dtype())
                        df_banemr['BAN_NB_F_FIN'] = round(df_banemr.BAN_NB_F_FIN * df_banemr.unit).astype(pd.Int64Dtype())
                        del df_banemr['BAN_FG_UNITE'], df_banemr['unit']
                        if coalesce:
                            df_banemr = coalesce_freqs(df_banemr) # A little bit long to compute, but this removes ~24000 useless/duplicate entries
                        #df.rename(columns={'BAN_ID': 'BANEMR_ID'}, inplace=True)
                        #df.set_index('BANEMR_ID', inplace=True)
                        df= df_banemr[['BAN_NB_F_DEB', 'BAN_NB_F_FIN']].drop_duplicates().sort_values(by='BAN_NB_F_DEB').reset_index(drop=True)
                        #df.to_sql('anfr_id_band', conn, if_exists='replace', index_label='BAN_ID', dtype={'BAN_ID': 'INTEGER primary key'})
                        df['BANSTR'] = df['BAN_NB_F_DEB'].astype(str) + '_' + df['BAN_NB_F_FIN'].astype(str)
                        df['BAN_ID'] = df.index
                        df.set_index('BANSTR', inplace=True)
                        df_banemr['BANSTR'] = df_banemr['BAN_NB_F_DEB'].astype(str) + '_' + df_banemr['BAN_NB_F_FIN'].astype(str)
                        df_banemr['BAN_ID'] = df.loc[df_banemr.BANSTR].BAN_ID.to_numpy()
                        del df_banemr['BANSTR'], df_banemr['BAN_NB_F_DEB'], df_banemr['BAN_NB_F_FIN']
                        df_banemr.to_sql('anfr_banemr', conn, if_exists='replace', index=False)
                        #df.reset_index(inplace=True)
                        df.set_index('BAN_ID', inplace=True)
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
    #create_views(dbfilename)

def gen_sites(dbfilename):
    with sqlite3.connect(dbfilename) as conn:
        cur = conn.cursor()
        global SYSLIST
        SYSLIST = [k[0] for k in cur.execute('select EMR_LB_SYSTEME from anfr_id_systeme').fetchall()]
        conn.create_function("mask_low", 1, mask_from_list_low64)
        conn.create_function("mask_high", 1, mask_from_list_high64)
        print("gen_sites")
        cur.executescript("""drop table if exists gen_sites ;
                             create table gen_sites as
                                select anfr_site.site_id, anfr_site.cor_dms, cor_nb_lon lon, cor_nb_lat lat,
                                count(distinct anfr_support.sup_id) support_count,group_concat(distinct anfr_support.sup_id) support_list, max(sup_nm_haut) h_max,
                                count(distinct sta_nm_anfr) sta_count, group_concat(distinct sta_nm_anfr) sta_list,
                                count(distinct anfr_antenne.aer_id) aer_count, group_concat(distinct anfr_antenne.aer_id) aer_list,
                                count(distinct emr_lb_systeme) tech_count, group_concat(distinct emr_lb_systeme) tech_list,
                                count(distinct anfr_emetteur.emr_id) emr_count, group_concat(distinct anfr_emetteur.emr_id) emr_list,
                                count(distinct ban_id) band_count,group_concat(distinct ban_id) band_list,
                                com_cd_insee insee, 0 as bitmask1, 0 as bitmask2

                                from anfr_site
                                inner join anfr_antenne on anfr_antenne.sup_id=anfr_support.sup_id
                                inner join anfr_support on anfr_site.SITE_ID=anfr_support.site_id
                                inner join anfr_emetteur on anfr_emetteur.aer_id=anfr_antenne.aer_id
                                inner join anfr_id_systeme on anfr_emetteur.sys_id=anfr_id_systeme.sys_id
                                inner join anfr_banemr on anfr_banemr.emr_id=anfr_emetteur.emr_id
                                group by anfr_site.site_id;

                             update gen_sites set bitmask1=mask_low(tech_list);
                             update gen_sites set bitmask2=mask_high(tech_list);
                          """) # anfr_site.cor_dms, nat_lb_nom, tpo_lb
                                #inner join anfr_id_proprietaire on anfr_id_proprietaire.tpo_id=anfr_support.tpo_id
                                #inner join anfr_id_nature on anfr_id_nature.nat_id=anfr_support.nat_id
                                #, group_concat(adm_lb_nom),
                                #inner join anfr_station on anfr_station.sta_nm_anfr=anfr_emetteur.sta_nm_anfr
                                #inner join anfr_id_exploitant on anfr_id_exploitant.adm_id=anfr_station.adm_id
                                #inner join anfr_id_nature on anfr_id_nature.nat_id=anfr_support.nat_id
                                #inner join anfr_id_proprietaire on anfr_id_proprietaire.tpo_id=anfr_support.tpo_id
                                # mask_low(group_concat(distinct emr_lb_systeme)), mask_high(group_concat(distinct emr_lb_systeme)) tech_bitmask2,
                                #inner join anfr_id_proprietaire on anfr_id_proprietaire.tpo_id=anfr_support.tpo_id

        print("gen_sectors")
        cur.executescript("""drop table if exists gen_sectors ;
                             create table gen_sectors as
                                select anfr_antenne.aer_id, aer_nb_azimut, aer_nb_alt_bas, cor_nb_lon lon, cor_nb_lat lat, anfr_antenne.sup_id,
                                sup_nm_haut h, anfr_emetteur.sta_nm_anfr,adm_lb_nom,
                                count(distinct emr_lb_systeme) tech_count, group_concat(distinct emr_lb_systeme) tech_list,
                                count(distinct anfr_emetteur.emr_id) emr_id_count, group_concat(distinct anfr_emetteur.emr_id) emr_id_list,
                                count(distinct ban_id) band_count,group_concat(distinct ban_id) band_list,
                                0 as bitmask1, 0 as bitmask2

                                from anfr_antenne
                                inner join anfr_support on anfr_antenne.sup_id=anfr_support.sup_id
                                inner join anfr_site on anfr_site.SITE_ID=anfr_support.site_id
                                inner join anfr_emetteur on anfr_emetteur.aer_id=anfr_antenne.aer_id
                                inner join anfr_id_systeme on anfr_emetteur.sys_id=anfr_id_systeme.sys_id
                                inner join anfr_banemr on anfr_banemr.emr_id=anfr_emetteur.emr_id
                                inner join anfr_station on anfr_station.sta_nm_anfr=anfr_emetteur.sta_nm_anfr
                                inner join anfr_id_exploitant on anfr_id_exploitant.adm_id=anfr_station.adm_id
                                group by anfr_antenne.aer_id;

                             update gen_sectors set bitmask1=mask_low(tech_list);
                             update gen_sectors set bitmask2=mask_high(tech_list);
                          """) # tae_lb,  tpo_lb,com_cd_insee insee,nat_lb_nom,
                                #inner join anfr_id_proprietaire on anfr_id_proprietaire.tpo_id=anfr_support.tpo_id
                                #inner join anfr_id_nature on anfr_id_nature.nat_id=anfr_support.nat_id
                                #inner join anfr_id_type_antenne on anfr_id_type_antenne.tae_id=anfr_antenne.tae_id
                                #mask_low(group_concat(distinct emr_lb_systeme)) tech_bitmask1, mask_high(group_concat(distinct emr_lb_systeme)) tech_bitmask2,

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
        #create_imtsectors_table(args.dbfile)
        #mask_v_support(args.dbfile)
        gen_sites(args.dbfile)


#################################################
# Legacy

IMT_MASK = ['GSM 900', 'GSM 1800', 'UMTS 900', 'UMTS 2100',
                'LTE 700', 'LTE 800', 'LTE 900', 'LTE 1800', 'LTE 2100', 'LTE 2600',
                '5G NR 700', '5G NR 1800', '5G NR 2100', '5G NR 3500']
IMT_OPS = ['ORANGE', 'SFR', 'FREE MOBILE', 'BOUYGUES TELECOM']

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
                        group by anfr_antenne.aer_id;


                        create table t_rop as
                        select anfr_site.site_id, anfr_antenne.sup_id, anfr_antenne.aer_id, aer_nb_azimut azh, aer_nb_alt_bas h, cor_dms, cor_nb_lat as lat, cor_nb_lon as lon, com_cd_insee insee , sta_nm_anfr_list stations, emr_id_list emetteurs,tech_list techs, op_list ops, sta_nm_anfr_count, emr_id_count, tech_count, op_count
                        from anfr_antenne
                        inner join (
                            select anfr_antenne.aer_id as aid,count(anfr_staant.sta_nm_anfr) sta_nm_anfr_count, group_concat(anfr_staant.sta_nm_anfr) sta_nm_anfr_list ,group_concat(adm_lb_nom) op_list, count(adm_lb_nom) op_count
                            from anfr_antenne inner join anfr_staant on anfr_staant.aer_id=anfr_antenne.aer_id
                            inner join anfr_station on anfr_station.sta_nm_anfr=anfr_staant.STA_NM_ANFR
                            inner join anfr_id_exploitant on anfr_station.adm_id=anfr_id_exploitant.adm_id
                            group by anfr_antenne.aer_id
                        ) on aid = anfr_antenne.aer_id
                        inner join (
                            select anfr_antenne.aer_id as aid2,count(emr_id) emr_id_count, group_concat(emr_id) emr_id_list,
                            count(distinct emr_lb_systeme) tech_count, group_concat(distinct emr_lb_systeme) tech_list
                            from anfr_antenne inner join anfr_emetteur on anfr_emetteur.aer_id=anfr_antenne.aer_id
                            inner join anfr_id_systeme on anfr_emetteur.sys_id=anfr_id_systeme.sys_id
                            group by anfr_antenne.aer_id
                        ) on aid2 = anfr_antenne.aer_id
                        inner join anfr_support on anfr_support.sup_id=anfr_antenne.sup_id
                        inner join anfr_site on anfr_site.SITE_ID=anfr_support.site_id
                        group by anfr_antenne.aer_id
                        having tech_list like 'GSM %' or tech_list like 'UMTS %' or tech_list like 'LTE %' or tech_list like '5G %'
                        """)

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
            (ban_nb_f_fin-ban_nb_f_deb)/1e3 bw_mhz,
            aer_nb_azimut azh, com_cd_insee insee, cor_dms as dms,
            cor_nb_lat as lat, cor_nb_lon as lon, emr_lb_systeme tech,
            adm_lb_nom op, tpo_lb op2, nat_lb_nom as sup_type
            from anfr_bande
            inner join anfr_id_band on anfr_bande.ban_id=anfr_id_band.BAN_ID
            inner join anfr_emetteur on anfr_bande.emr_id=anfr_emetteur.emr_id
            inner join anfr_antenne on anfr_emetteur.aer_id=anfr_antenne.aer_id
            inner join anfr_support on anfr_antenne.sup_id=anfr_support.sup_id
            inner join anfr_station on anfr_antenne.sta_nm_anfr=anfr_station.sta_nm_anfr
            inner join anfr_id_exploitant on anfr_station.adm_id=anfr_id_exploitant.adm_id
            inner join anfr_id_proprietaire on anfr_support.tpo_id=anfr_id_proprietaire.tpo_id
            inner join anfr_id_nature on anfr_support.nat_id=anfr_id_nature.nat_id
            inner join anfr_id_systeme on anfr_emetteur.sys_id=anfr_id_systeme.sys_id
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


def mask_v_support_legacy(dbfilename):
    with sqlite3.connect(dbfilename) as conn:
        cur = conn.cursor()
        global SYSLIST
        SYSLIST = [k[0] for k in cur.execute('select EMR_LB_SYSTEME from anfr_id_systeme').fetchall()]
        conn.create_function("mask_low", 1, mask_from_list_low64)
        conn.create_function("mask_high", 1, mask_from_list_high64)
        # FIXME: better create gen_sites ?
        cur.executescript("""drop table if exists gen_support ;
                             create table gen_support as
                                    select anfr_support.sup_id, anfr_support.site_id, cor_nb_lat as lat, cor_nb_lon as lon, sup_nm_haut, nat_lb_nom, tpo_lb, sta_nm_anfr_count,sta_nm_anfr_list, aer_id_count,aer_id_list,com_cd_insee insee,mask_low(tech_list) tech_bitmask1, mask_high(tech_list) tech_bitmask2, tech_count,tech_list,cor_dms
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
                                    group by anfr_support.sup_id;""")

#df = coalesce_freqs(df.groupby("EMR_ID"))
# def coalesce_freqs_slow(dfg):
#     res = []
#     kk=0
#     for EMR_ID,freqs in dfg:
#         sys.stderr.write(f"\rCoalescing freq entries: {100 * kk // len(dfg)} %")
#         kk+=1
#         groups = []
#         for row in freqs.iterrows():
#             fmin = row[1].BAN_NB_F_DEB
#             fmax = row[1].BAN_NB_F_FIN
#             if pd.isna(fmin) or pd.isna(fmax):
#                 continue
#             found=False
#             for g in groups:
#                 gmin = g[0] ; gmax = g[1]
#                 if fmin>=gmin and fmin<=gmax and fmax>=gmin and fmax<=gmax: found=True; break # Already covered
#                 elif fmin>=gmin and fmin<=gmax: g[1] = fmax ; found=True ; break # Extend upwards
#                 elif fmax>=gmin and fmax<=gmax: g[0] = fmin ; found=True ; break # Extend downwards
#             if found==False: groups.append([fmin,fmax,EMR_ID])
#         res.extend(groups)
#     return pd.DataFrame(data=res, columns=["fmin","fmax","EMR_ID"])

