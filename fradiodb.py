#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Adrien DEMAREZ
"""

import sqlite3
from glob import glob
from os import sep,makedirs,remove,getcwd
from os.path import splitext,basename,exists
import zipfile
import csv
import io
import argparse
import urllib.request as req
import struct
import json
from http.server import HTTPServer, BaseHTTPRequestHandler
import re

# Schema diagram: d2 ./schema.d2 -s -t 100

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

def wkb_point(lon, lat, srs_id=4326):
    # https://www.geopackage.org/spec/#gpb_data_blob_format
    magic = b"GP"            # magic
    version = 0              # version
    flags = 1                # little endian, no envelope
    header = struct.pack("<2sBBi", magic, version, flags, srs_id)

    # WKB geometry
    byte_order = 1           # little endian
    geom_type = 1            # POINT
    wkb = struct.pack("<BIdd", byte_order, geom_type, lon, lat)
    return header + wkb

def import_anfr_zip(dbfilename, dirpath='anfr'):
    if exists(dbfilename):
        print("Deleting previous DB")
        remove (dbfilename)
    with sqlite3.connect(dbfilename) as conn:
        conn.create_function("pytitle", 1, str.title)
        conn.create_function("pyfix_cp", 2, lambda x, y: -1)
        conn.create_function("pyfix_insee", 2, lambda x, y: -1)
        conn.create_function("pyprint", 1, lambda x: print(x))
        conn.create_function("pycoalesce", 1, coalesce_freqs)
        conn.create_function("pywbkpoint", 2, wkb_point)
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

def get_syslist(dbfilename):
    with sqlite3.connect(dbfilename) as conn:
        cur = conn.cursor()
        SYSLIST = [k[0] for k in cur.execute('select system from id_systems').fetchall()]
    return SYSLIST

def mask_from_list(strlist, masklist=None):
    if strlist is None: return None
    if masklist is None: masklist = SYSLIST
    mask = 0
    for entry in strlist.split(','):
        mask |= (1<<masklist.index(entry))
    return mask
# MAX_BIT=52+1 because 1°/ sqlite does not deal well with unsigned 64 ints so the last bit is not usable, and
#   2°/ Excel and Libreoffice have accuracy issues with MAX_BIT>52 because they compute as IEEE754 double floats and the mantissa has 52 useful bits (and I want to keep the spreadsheet for easy mask calculations)
MAX_BIT = 52+1
def mask128_low64(mask128): return mask128 & ((1<<MAX_BIT)-1)
def mask128_high64(mask128): return mask128 >> MAX_BIT
def mask_from_list_low64(strlist, masklist=None): return mask128_low64(mask_from_list(strlist, masklist))
def mask_from_list_high64(strlist, masklist=None): return mask128_high64(mask_from_list(strlist, masklist))

def masks64_to_mask128(mask_low, mask_high): return mask_low | (mask_high << MAX_BIT)
def list_from_mask(mask, masklist=None):
    if masklist is None: masklist = SYSLIST
    alist = []
    for bit in range(len(masklist)):
        if mask & (1<<bit):
            alist.append(masklist[bit])
    return ','.join(alist)

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

def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d

def query(dbfilename, itemid, itemtype="site"):
    conn = sqlite3.connect(dbfilename)
    conn.row_factory = dict_factory
    if itemtype=="site":
        SYSLIST = get_syslist(dbfilename) # FIXME: cache this result
        return json.dumps(query_site(conn, itemid, SYSLIST))
    elif itemtype=='support':
        return json.dumps(query_antennas(conn, itemid))
    elif itemtype=='antenna':
        return json.dumps(query_transmitters(conn, itemid))

def query_site(conn, site_id, syslist):
    cur  = conn.cursor()
    res = cur.execute("""select dms, lon, lat, address, postcode, inseecode, tech_bitmask1, tech_bitmask2 from sites where sites.id=?""", (site_id,)).fetchall()[0]
    # res['stations'] = cur.execute("""
    #     select cast(group_concat(distinct station_id) as text) sta_list,
    #     from sites
    #     inner join antennas on antennas.sup_id=supports.id
    #     inner join supports on sites.id=supports.site_id
    #     inner join transmitters on transmitters.antenna_id=antennas.id
    #     inner join stations on stations.id=transmitters.station_id where sites.id=?
    #     group by sites.id""", (site_id,)).fetchall()[0]
    res['tech_list'] = list_from_mask(masks64_to_mask128(res['tech_bitmask1'],res['tech_bitmask2']),syslist)
    del res['tech_bitmask1'],res['tech_bitmask2']
    res['supports'] = query_supports(conn, site_id)
    return res

def query_supports(conn, site_id):
    cur  = conn.cursor()
    res_supports = cur.execute("""select supports.id support_id, support_type, support_owner, sup_height support_height from supports
                                         inner join id_support_owners on id_support_owners.id=owner_id
                                         inner join id_support_types on id_support_types.id=suptype_id
                                         where site_id=?""", (site_id,)).fetchall()
    for res in res_supports:
        res['antennas'] = query_antennas(conn, res['support_id'])
    return res_supports

def query_antennas(conn, support_id):
    cur  = conn.cursor()
    res_antennas = cur.execute("""select antennas.id antenna_id, dimension, azimuth, ant_height, dim_type, antenna_type from antennas
                                         inner join id_antenna_types on anttype_id=id_antenna_types.id
                                         where sup_id=?""", (support_id,)).fetchall()
    for res in res_antennas:
        res['transmitters'] = query_transmitters(conn, res['antenna_id'])
    return res_antennas

def query_transmitters(conn, antenna_id):
    cur  = conn.cursor()
    res_transmitters = cur.execute("""select station_name, operator, system, date(date_switchedon,'unixepoch') date_on, bandstr band_kHz
                                         from transmitters
                                         inner join id_systems on system_id=id_systems.id
                                         inner join stations on station_id=stations.id
                                         inner join id_operators on operator_id=id_operators.id
                                         inner join v_bandgroups on bandgroup_id=v_bandgroups.id
                                         where antenna_id=?""", (antenna_id,)).fetchall()
    return res_transmitters 

class fradiodb_web_handler(BaseHTTPRequestHandler):
    def staticfile(self, mime):
        enc = self.headers.get("Accept-Encoding", "")
        encoding = None
        if "gzip" in enc and exists(getcwd()+self.path+'.gz'):
            self.path += ".gz"
            encoding = "gzip"
        try:
            with open(getcwd()+self.path, "rb") as f:
                content = f.read()
            self.send_response(200)
            self.send_header("Content-Type", mime)
            if encoding:
                self.send_header("Content-Encoding", encoding)
            self.end_headers()
            self.wfile.write(content)
        except FileNotFoundError:
            self.send_error(404)
        return

    def do_GET(self):
        if self.path=="/": self.path="/index.html"
        if self.path == "/index.html":
            return self.staticfile("text/html")
        if self.path == "/minmax.json":
            return self.staticfile("application/json")
        if self.path == "/anfr.json":
            return self.staticfile("application/geo+json")
        if self.path == "/anfr.fgb":
            return self.staticfile("application/x-flatgeobuf")

        # Dynamic route /site/<site_id>
        match = re.match(r"^/site/([^/]+)$", self.path)
        if match:
            site_id = match.group(1)
            payload = query(self.server.dbfile, int(site_id), "site")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.end_headers()
            self.wfile.write(payload.encode("utf-8"))
            return

        # Otherwise 404
        self.send_error(404)

def gen_minmax_json(dbfilename="anfr2026.gpkg", jsonfile="minmax.json"):
    conn = sqlite3.connect(dbfilename)
    conn.row_factory = dict_factory
    query = """select min(support_count) sup_c_min, max(support_count) sup_c_max, min(h) h_min, max(h) h_max,
               min(operator_count) op_c_min, max(operator_count) op_c_max, min(sta_count) sta_c_min, max(sta_count) sta_c_max,
               min(ant_count) ant_c_min, max(ant_count) ant_c_max
                   from (select count(distinct supports.id) support_count, max(sup_height) h, count(distinct operator) operator_count,
                    count(distinct station_id) sta_count,count(distinct antennas.id) ant_count
                    from sites
                    inner join antennas on antennas.sup_id=supports.id
                    inner join supports on sites.id=supports.site_id
                    inner join transmitters on transmitters.antenna_id=antennas.id
                    inner join id_systems on transmitters.system_id=id_systems.id
                    inner join stations on stations.id=transmitters.station_id
                    inner join id_operators on id_operators.id=stations.operator_id
                    group by sites.id)"""
    cur = conn.cursor()
    res = cur.execute(query).fetchall()[0]
    metrics = {
        "h":  [res['h_min'], res['h_max']],
        #"support_count": [res['sup_c_min'], res['sup_c_max']],
        #"operator_count": [res['op_c_min'], res['op_c_max']],
        "sta_count": [res['sta_c_min'], res['sta_c_max']],
        "ant_count": [res['ant_c_min'], res['ant_c_max']],
    }
    with open(jsonfile, "w") as f:
        json.dump(metrics, f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dbfile", help="DB path")
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    parser_create = subparsers.add_parser('create', help="Create DB")
    parser_create.add_argument("--download", "-d", help="Download data from data.gouv.fr", action='store_true', default=False)
    parser_create.add_argument("--datadir", "-p", help="Data location", default="anfr")

    parser_query = subparsers.add_parser('info', help="Query DB")
    parser_query.add_argument("id", help="ID to query")
    parser_query.add_argument("--type", "-t", help="site | support | station | antenna | transmitter", default="site")

    parser_server = subparsers.add_parser('serve', help="REST server")
    parser_server.add_argument("--port", "-p", help="Port", default=8001)

    args = parser.parse_args()

    if args.subcommand=='create':
        if args.download:
            download_data(args.datadir)
    #if args.importdb is not None:
        #import_cities(args.dbfile, dirpath=mydir)
        import_anfr_zip(args.dbfile, dirpath=args.datadir)
    elif args.subcommand=='info':
        res = query(args.dbfile, args.id, args.type)
        print(res)
    elif args.subcommand=='serve':
        server = HTTPServer(("localhost", int(args.port)), fradiodb_web_handler)
        server.dbfile = args.dbfile
        print(f"Server running on http://localhost:{args.port}")
        server.serve_forever()
