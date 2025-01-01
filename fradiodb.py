#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Adrien DEMAREZ
"""

import sqlite3
from glob import glob
from os import sep,makedirs,remove
from os.path import splitext,basename,exists
import zipfile
import csv
import io
import argparse
import urllib.request as req

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

def coalesce_freqs(bandlist_str):
    # bandlist_str is a string listing several bands fmin_fmax, sorted by fmin
    # small trick: (fmin,fmax) is encoded as (real,imag) to avoid nested lists and associated performance issues
    if bandlist_str is None: return None
    if not ',' in bandlist_str: return bandlist_str
    bands = [0+0j]
    for band in bandlist_str.replace('[','').replace(']','').replace('"','').split(','):
        fminmax = band.split('_')
        fmin = int(float(fminmax[0]))
        fmax = int(float(fminmax[1]))
        gmin = int(bands[-1].real)
        gmax = int(bands[-1].imag)
        if fmin>=gmin and fmin<=gmax and fmax>=gmin and fmax<=gmax: pass  # Band included within previous one
        elif fmin>=gmin and fmin<=gmax: bands[-1] = gmin+1j*fmax  # Extend previous band upwards
        elif fmax>=gmin and fmax<=gmax: # Extend previous band downwards. FIXME: should not happen if list sorted by fmin
            bands[-1] = fmin+1j*gmax
            print(f"error {bands}")
        else: bands.append(fmin + 1j*fmax)
    res = '["' + '","'.join([str(band.real) + '_' + str(band.imag) for band in bands[1:]]) + '"]'
    return res 

def import_anfr_zip(dbfilename, dirpath='anfr'):
    if exists(dbfilename):
        remove (dbfilename)
    with sqlite3.connect(dbfilename) as conn:
        conn.create_function("pytitle", 1, str.title)
        conn.create_function("pyfix_cp", 2, lambda x, y: -1)
        conn.create_function("pyfix_insee", 2, lambda x, y: -1)
        conn.create_function("pyprint", 1, lambda x: print(x))
        conn.create_function("pycoalesce", 1, coalesce_freqs)
        cur = conn.cursor()

        # Import CSVs into SQL tables
        for myzipfile in glob(dirpath + sep + "*anfr*.zip"):
            with zipfile.ZipFile(myzipfile) as zFile:
                for csvfile in zFile.infolist():
                    tablename = splitext(basename(csvfile.filename))[0].lower()
                    with zFile.open(csvfile) as binz, io.TextIOWrapper(binz, encoding='iso8859-1' if tablename == "sup_support" else "utf-8") as textz:
                        print("importing " + csvfile.filename)
                        reader = csv.reader(textz, delimiter=";")
                        headers = next(reader)
                        cur.execute(f"create temp table {tablename} ({','.join(headers)})")
                        cur.executemany(f"insert into {tablename} values ({','.join(['?' for _ in headers])});", reader)
                        for h in headers:
                            cur.execute(f"update {tablename} set {h}=null where trim({h}) in ('','-',char(9))")

        # Perform the main transformations (see SQL script)
        with open("schema.sql", encoding='utf-8') as schema:
            cur.executescript(schema.read())

        # Compute tech_bitmask
        global SYSLIST
        SYSLIST = [k[0] for k in cur.execute('select system from id_systems').fetchall()]
        conn.create_function("mask_low", 1, mask_from_list_low64)
        conn.create_function("mask_high", 1, mask_from_list_high64)
        for table in ('sites', 'supports', 'stations', 'antennas'):
            print(f'Computing tech_bitmask for {table}')
            cur.executescript(f"""
                create temp table tmp_{table} as
                    select {table}.id tmp_id, group_concat(distinct system) tech_list
                    from {table}
                    {"inner join supports on sites.id=supports.site_id" if table=="sites" else ""}
                    {"inner join antennas on antennas.sup_id=supports.id" if table in ("sites", "supports") else ""}
                    inner join transmitters on transmitters.{"station_id=stations.id" if table=="stations" else "antenna_id=antennas.id"}
                    inner join id_systems on transmitters.system_id=id_systems.id
                    group by {table}.id;
                update {table} set tech_bitmask1=mask_low(tech_list) from (select tmp_id, tech_list from tmp_{table}) foo where foo.tmp_id={table}.id;
                update {table} set tech_bitmask2=mask_high(tech_list) from (select tmp_id, tech_list from tmp_{table}) foo where foo.tmp_id={table}.id;
                drop table tmp_{table};
            """)

def mask_from_list_low64(strlist, masklist=None):
    if strlist is None: return None
    if masklist is None: masklist = SYSLIST
    mask = 0
    for entry in strlist.split(','):
        mask |= (1<<masklist.index(entry))
    return mask128_low64(mask)
def mask_from_list_high64(strlist, masklist=None):
    if strlist is None: return None
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

    args = parser.parse_args()
    if args.download:
        download_data(args.datadir)
    if args.importdb is not None:
        #import_cities(args.dbfile, dirpath=mydir)
        import_anfr_zip(args.importdb, dirpath=args.datadir)
