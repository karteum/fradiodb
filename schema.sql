-- MIME: application/geopackage+sqlite3
-- Required GeoPackage SQLite settings
PRAGMA application_id = 1196437808;   -- 0x47504B47 ("GPKG")
PRAGMA user_version = 10300;          -- GeoPackage version 1.3.0

-- Spatial reference systems
CREATE TABLE gpkg_spatial_ref_sys (
    srs_name TEXT NOT NULL,
    srs_id INTEGER PRIMARY KEY,
    organization TEXT NOT NULL,
    organization_coordsys_id INTEGER NOT NULL,
    definition  TEXT NOT NULL,
    description TEXT
);
INSERT INTO gpkg_spatial_ref_sys VALUES ('Undefined Cartesian SRS', -1, 'NONE', -1, 'undefined', 'undefined Cartesian coordinate reference system');
INSERT INTO gpkg_spatial_ref_sys VALUES ('Undefined geographic SRS', 0, 'NONE', 0, 'undefined', 'undefined geographic coordinate reference system');
INSERT INTO gpkg_spatial_ref_sys VALUES ('WGS 84 geodetic', 4326, 'EPSG', 4326, 'GEOGCS["WGS 84",DATUM["World Geodetic System 1984", SPHEROID["WGS 84",6378137,298.257223563]], PRIMEM["Greenwich",0], UNIT["degree",0.0174532925199433]]', 'longitude/latitude coordinates in decimal degrees on the WGS 84 spheroid');

-- Contents table
CREATE TABLE gpkg_contents (
    table_name TEXT NOT NULL PRIMARY KEY,
    data_type TEXT NOT NULL,
    identifier TEXT UNIQUE,
    description TEXT DEFAULT '',
    last_change DATETIME NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%fZ','now')),
    min_x DOUBLE,
    min_y DOUBLE,
    max_x DOUBLE,
    max_y DOUBLE,
    srs_id INTEGER,
    CONSTRAINT fk_gc_r_srs_id FOREIGN KEY (srs_id) REFERENCES gpkg_spatial_ref_sys(srs_id)
);
INSERT INTO gpkg_contents (table_name, data_type, identifier, description, srs_id) VALUES ('sites', 'features', 'sites', 'Radio tower locations', 4326);
INSERT INTO gpkg_contents (table_name, data_type, identifier) VALUES ('stations', 'attributes', 'stations');
INSERT INTO gpkg_contents (table_name, data_type, identifier) VALUES ('antennas', 'attributes', 'antennas');
INSERT INTO gpkg_contents (table_name, data_type, identifier) VALUES ('bandgroups', 'attributes', 'bandgroups');
INSERT INTO gpkg_contents (table_name, data_type, identifier) VALUES ('transmitters', 'attributes', 'transmitters');
INSERT INTO gpkg_contents (table_name, data_type, identifier) VALUES ('supports', 'attributes', 'supports');
INSERT INTO gpkg_contents (table_name, data_type, identifier) VALUES ('id_systems', 'attributes', 'id_systems');
INSERT INTO gpkg_contents (table_name, data_type, identifier) VALUES ('id_antenna_types', 'attributes', 'id_antenna_types');
INSERT INTO gpkg_contents (table_name, data_type, identifier) VALUES ('id_operators', 'attributes', 'id_operators');
INSERT INTO gpkg_contents (table_name, data_type, identifier) VALUES ('id_support_owners', 'attributes', 'id_support_owners');
INSERT INTO gpkg_contents (table_name, data_type, identifier) VALUES ('id_support_types', 'attributes', 'id_support_types');

-- Geometry columns table
CREATE TABLE gpkg_geometry_columns (
    table_name TEXT NOT NULL,
    column_name TEXT NOT NULL,
    geometry_type_name TEXT NOT NULL,
    srs_id INTEGER NOT NULL,
    z TINYINT NOT NULL,
    m TINYINT NOT NULL,
    CONSTRAINT pk_geom_cols PRIMARY KEY (table_name, column_name),
    CONSTRAINT uk_gc_table_name UNIQUE (table_name),
    CONSTRAINT fk_gc_tn FOREIGN KEY (table_name) REFERENCES gpkg_contents(table_name),
    CONSTRAINT fk_gc_srs FOREIGN KEY (srs_id) REFERENCES gpkg_spatial_ref_sys (srs_id)
);
INSERT INTO gpkg_geometry_columns (table_name, column_name, geometry_type_name, srs_id, z, m) VALUES ('sites', 'geom', 'POINT', 4326, 0, 0);


CREATE TABLE sample_feature_table (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  geometry GEOMETRY,
  text_attribute TEXT,
  real_attribute REAL,
  boolean_attribute BOOLEAN,
  raster_or_photo BLOB
);
CREATE TABLE gpkg_tile_matrix_set (
  table_name TEXT NOT NULL PRIMARY KEY,
  srs_id INTEGER NOT NULL,
  min_x DOUBLE NOT NULL,
  min_y DOUBLE NOT NULL,
  max_x DOUBLE NOT NULL,
  max_y DOUBLE NOT NULL,
  CONSTRAINT fk_gtms_table_name FOREIGN KEY (table_name) REFERENCES gpkg_contents(table_name),
  CONSTRAINT fk_gtms_srs FOREIGN KEY (srs_id) REFERENCES gpkg_spatial_ref_sys (srs_id)
);
CREATE TABLE gpkg_tile_matrix (
  table_name TEXT NOT NULL,
  zoom_level INTEGER NOT NULL,
  matrix_width INTEGER NOT NULL,
  matrix_height INTEGER NOT NULL,
  tile_width INTEGER NOT NULL,
  tile_height INTEGER NOT NULL,
  pixel_x_size DOUBLE NOT NULL,
  pixel_y_size DOUBLE NOT NULL,
  CONSTRAINT pk_ttm PRIMARY KEY (table_name, zoom_level),
  CONSTRAINT fk_tmm_table_name FOREIGN KEY (table_name) REFERENCES gpkg_contents(table_name)
);

------------------------------------------------------------------------------------------------
-- Starting point : tables are imported from CSV (table name as filename, all fields as text)
--   N.B. the following script uses user-defined functions in the calling Python code (i.e. in case the calling code is not the original Python script, those functions need to be implemented or the code below has to be modified accordingly):
--   * pytitle (Python's str.title(), may be replaced by lower())
--   * pyfix_cp/pyfix_insee (correction of wrong zipcodes, see below)
--   * pyprint (convenient for debugging, easy to remove if needed)
--   * pycoalesce (regroup/coalesce entries for adjacent frequency bands)

select pyprint('Processing ID tables');
create table id_antenna_types (id integer primary key, antenna_type text) strict;
insert into id_antenna_types select TAE_ID, TAE_LB from sup_type_antenne;

create table id_support_types (id integer primary key, support_type text) strict;
insert into id_support_types select NAT_ID, NAT_LB_NOM from sup_nature;

create table id_support_owners (id integer primary key, support_owner text) strict;
insert into id_support_owners select TPO_ID, TPO_LB from sup_proprietaire;

create table id_operators (id integer primary key, operator text) strict;
insert into id_operators select ADM_ID, ADM_LB_NOM from sup_exploitant;

select pyprint('Processing sup_station');
update sup_station set DTE_EN_SERVICE=NULL where DTE_EN_SERVICE='';
update sup_station set DTE_MODIF=NULL where DTE_MODIF='';
create table stations (
    id integer primary key,
    operator_id integer,
    station_name text unique,
    date_created integer,
    date_modified integer,
    date_operational integer,
    tech_bitmask1 integer,
    tech_bitmask2 integer,
    foreign key(operator_id) references id_operators(id)
) strict;
insert into stations select DEM_NM_COMSIS, ADM_ID, STA_NM_ANFR, NULL, NULL, -- store dates as unix epoch
    strftime('%s', substr(DTE_IMPLANTATION, 7, 4) || '-' || substr(DTE_IMPLANTATION, 4, 2) || '-' || substr(DTE_IMPLANTATION, 1, 2)),
    strftime('%s', substr(DTE_MODIF, 7, 4) || '-' || substr(DTE_MODIF, 4, 2) || '-' || substr(DTE_MODIF, 1, 2)),
    strftime('%s', substr(DTE_EN_SERVICE, 7, 4) || '-' || substr(DTE_EN_SERVICE, 4, 2) || '-' || substr(DTE_EN_SERVICE, 1, 2))
    from sup_station;

select pyprint('Processing sup_support');
-- Several transformations made here in order to
--   1°/ have latitude and longitude as floats (but we keep a "dms" string which is useful as a key for further processing)
--   2°/ regroup the 4 address fields in a single one, together with some cleanups
--   3°/ regroup all data w.r.t. a site (i.e. address, lat, lon, zipcode, etc) into a dedicated table "sites" instead of repeating it for multiple supports (several supports may be at the same address). Notice that in some cases, the same dms location was leading to various zipcodes (this is corrected with pyfix_cp() and pyfix_insee()) or to different addresses (sometimes similar and sometimes not : in those cases the distinct strings are aggregated with ' ¤ ')
--   4°/ also avoid the repetitions in table "support" (since one support may host several stations and one station may include equipment on several supports, the field "sta_nm_anfr" was leading to many repeated lines in this table)
-- select dms, group_concat(distinct adr_nm_cp), group_concat(distinct com_cd_insee), count(dms) c from (select distinct dms, adr_nm_cp, com_cd_insee from tmp_supports) foo group by dms having c>1
create temp table tmp_supports as select distinct
    sup_id, NAT_ID, SUP_NM_HAUT, tpo_id, ADR_NM_CP, COM_CD_INSEE, --sta_nm_anfr
    COR_NB_DG_LAT||'°'||COR_NB_MN_LAT||''''||COR_NB_SC_LAT||'"'||COR_CD_NS_LAT||' '||COR_NB_DG_LON||'°'||COR_NB_MN_LON||''''||COR_NB_SC_LON||'"'||COR_CD_EW_LON dms,
    round((cast(COR_NB_DG_LON as real) + cast(COR_NB_MN_LON as real) / 60 + cast(COR_NB_SC_LON as real) / 3600) * (case COR_CD_EW_LON when 'E' then 1 else -1 end), 4) lon,
    round((cast(COR_NB_DG_LAT as real) + cast(COR_NB_MN_LAT as real) / 60 + cast(COR_NB_SC_LAT as real) / 3600) * (case COR_CD_NS_LAT when 'N' then 1 else -1 end), 4) lat,
    pytitle(concat_ws(', ', trim(ADR_LB_LIEU), trim(ADR_LB_ADD1), trim(ADR_LB_ADD2), trim(ADR_LB_ADD3))) address
    from sup_support;
update tmp_supports set ADR_NM_CP=pyfix_cp(lon, lat), COM_CD_INSEE=pyfix_insee(lon, lat) where dms in
    (select dms from (select distinct dms, adr_nm_cp, com_cd_insee from tmp_supports)
    group by dms having count(dms)>1);

create table sites (
    id integer primary key,
    dms text unique,
    lon real,
    lat real,
    address text,
    postcode integer,
    inseecode text,
    tech_bitmask1 integer,
    tech_bitmask2 integer,
    geom BLOB
) strict;
insert into sites(dms,lon,lat,postcode,inseecode, geom) select distinct dms, lon, lat, ADR_NM_CP, COM_CD_INSEE, pywbkpoint(lon,lat) from tmp_supports;
update sites set address=addgroup from -- Unfortunately sqlite does not support group_concat(distinct address, ' ¤ ') so this is a workaround
    (select dms, group_concat(address, ' ¤ ') addgroup from
        (select distinct dms, address from tmp_supports)
    group by dms) foo -- having count(dms)>1
where sites.dms=foo.dms;

create table supports (
    id integer primary key,
    suptype_id integer,
    owner_id integer,
    site_id integer,
    sup_height integer,
    tech_bitmask1 integer,
    tech_bitmask2 integer,
    foreign key(suptype_id) references id_support_types(id),
    foreign key(owner_id) references id_support_owners(id),
    foreign key(site_id) references sites(id)
) strict;
insert into supports select distinct cast(sup_id as integer), cast(NAT_ID as integer), cast(TPO_ID as integer), sites.id, cast(SUP_NM_HAUT as real), NULL, NULL
    from tmp_supports inner join sites on tmp_supports.dms=sites.dms;
drop table tmp_supports;

select pyprint('Processing sup_antenna');
create table antennas (
    id integer primary key,
    anttype_id integer,
    sup_id integer,
    dimension real,
    azimuth integer,
    ant_height real,
    dim_type text,
    tech_bitmask1 integer,
    tech_bitmask2 integer,
    foreign key(anttype_id) references id_antenna_types(id),
    foreign key(sup_id) references supports(id)
) strict;
insert into antennas select distinct AER_ID, TAE_ID, SUP_ID, cast(AER_NB_DIMENSION as real), cast(AER_NB_AZIMUT as integer),
    cast(AER_NB_ALT_BAS as real), AER_FG_RAYON, NULL, NULL from sup_antenne; -- tech_bitmask shall be added later

select pyprint('Processing sup_bande');
-- Several transformations made here in order to factorize the data from sup_bande in a much shorter table "bandgroups", noting that one transmitter is associated with a single bandgroup and many transmitters operate with the same set of bands (e.g. sites from mobile operators)
create temp view tmp_bands as select distinct EMR_ID, round(fmin*unit) fmin_kHz, round(fmax*unit) fmax_kHz, cast(round(fmin*unit) as text) ||'_'||cast(round(fmax*unit) as text) fstr from
        (select EMR_ID, (case BAN_FG_UNITE when 'K' then 1 when 'M' then 1000 when 'G' then 1000000 end) unit, cast(BAN_NB_F_DEB as real) fmin, cast(BAN_NB_F_FIN as real) fmax
        from sup_bande where BAN_FG_UNITE is not null and BAN_NB_F_DEB is not null and BAN_NB_F_FIN is not null) -- FIXME: and BAN_NB_F_DEB!=BAN_NB_F_FIN ?
    order by fmin_kHz, fmax_kHz, EMR_ID;
create temp view tmp_bands2 as select emr_id, pycoalesce('["'||group_concat(fstr, '","')||'"]') jstr from tmp_bands group by emr_id;
create temp table tmp_bands3 as select distinct jstr from tmp_bands2; -- must be a table (not view) because we further need rowid
create temp view tmp_bands4 as select t.rowid id, j.value fstr from tmp_bands3 t, json_each(t.jstr) j; -- explode()
create temp view tmp_bands5 as select emr_id, tmp_bands3.rowid bandgroup_id from tmp_bands2 inner join tmp_bands3 on tmp_bands2.jstr=tmp_bands3.jstr;
create table bandgroups (
    id integer,
    fmin_kHz integer,
    fmax_kHz integer
) strict;
insert into bandgroups select id, cast(substr(fstr, 0, instr(fstr, '_')) as integer) fmin_kHz, cast(substr(fstr, 1+instr(fstr, '_')) as integer) fmax_kHz from tmp_bands4;

select pyprint('Processing sup_emetteur');
-- Group "system" label in a dedicated table (which both avoids repetitions, and enables to further compute tech_bitmask)
create table id_systems (id integer primary key, system text) strict;
insert into id_systems(system) select distinct EMR_LB_SYSTEME from sup_emetteur where EMR_LB_SYSTEME != '' order by EMR_LB_SYSTEME;

-- Notice that several transmitters may be associated to a single antenna (e.g. MORAN) and several stations may be associated to a single transmitter (e.g. MOCN)
create table transmitters (
    id integer primary key,
    station_id integer,
    antenna_id integer,
    system_id integer,
    bandgroup_id integer,
    date_switchedon integer,
    foreign key(station_id) references stations(id),
    foreign key(antenna_id) references antennas(id),
    foreign key(system_id) references systems(id),
    foreign key(bandgroup_id) references bandgroups(id)
) strict;
insert into transmitters select sup_emetteur.EMR_ID, stations.id, AER_ID, id_systems.id, bandgroup_id,
    strftime('%s', substr(EMR_DT_SERVICE, 7, 4) || '-' || substr(EMR_DT_SERVICE, 4, 2) || '-' || substr(EMR_DT_SERVICE, 1, 2))
    from sup_emetteur inner join id_systems on sup_emetteur.EMR_LB_SYSTEME=id_systems.system
    inner join stations on sup_emetteur.STA_NM_ANFR=stations.station_name
    inner join tmp_bands5 on tmp_bands5.emr_id=sup_emetteur.emr_id;

select pyprint('Adding views');
-- Add a couple of views for convenience
create view v_bandgroups as select id, count(fmin_kHz) band_count, group_concat(fmin_kHz||'-'||fmax_kHz) bandstr, sum(fmax_kHz-fmin_kHz) BWtot from bandgroups group by id;

create view v_sites as
    select sites.id, sites.dms, lon, lat, inseecode,
    count(distinct operator) operator_count, cast(group_concat(distinct operator) as text) operator_list,
    count(distinct supports.id) support_count, group_concat(distinct supports.id) support_list, max(sup_height) h_max,
    count(distinct station_id) sta_count, cast(group_concat(distinct station_id) as text) sta_list,
    count(distinct antennas.id) ant_count, cast(group_concat(distinct antennas.id) as text) ant_list,
    count(distinct system) tech_count, cast(group_concat(distinct system) as text) tech_list,
    count(distinct transmitters.id) tx_count, cast(group_concat(distinct transmitters.id) as text) tx_list,
    count(distinct bandstr) band_count, cast(group_concat(distinct bandstr) as text) band_list,
    sites.tech_bitmask1, sites.tech_bitmask2
    from sites
    inner join antennas on antennas.sup_id=supports.id
    inner join supports on sites.id=supports.site_id
    inner join transmitters on transmitters.antenna_id=antennas.id
    inner join id_systems on transmitters.system_id=id_systems.id
    inner join v_bandgroups on transmitters.bandgroup_id=v_bandgroups.id
    inner join stations on stations.id=transmitters.station_id
    inner join id_operators on id_operators.id=stations.operator_id
    group by sites.id;

create view v_antennas as
    select antennas.id, azimuth, ant_height, lon, lat, antennas.sup_id,
    sup_height, transmitters.station_id,
    count(distinct operator) operator_count, cast(group_concat(distinct operator) as text) operator_list,
    count(distinct system) tech_count, cast(group_concat(distinct system) as text) tech_list,
    count(distinct transmitters.id) tx_count, cast(group_concat(distinct transmitters.id) as text) tx_list,
    count(distinct bandstr) band_count, cast(group_concat(distinct bandstr) as text) band_list,
    antennas.tech_bitmask1, antennas.tech_bitmask2
    from antennas
    inner join supports on antennas.sup_id=supports.id
    inner join sites on sites.id=supports.site_id
    inner join transmitters on transmitters.antenna_id=antennas.id
    inner join id_systems on transmitters.system_id=id_systems.id
    inner join stations on stations.id=transmitters.station_id
    inner join v_bandgroups on transmitters.bandgroup_id=v_bandgroups.id
    inner join id_operators on id_operators.id=stations.operator_id
    group by antennas.id;

create view v_transmitters as
    select transmitters.id, antenna_id, station_name, operator, bandstr, system
    from transmitters
    inner join stations on stations.id=transmitters.station_id
    inner join id_operators on id_operators.id=stations.operator_id
    inner join id_systems on transmitters.system_id=id_systems.id
    inner join v_bandgroups on transmitters.bandgroup_id=v_bandgroups.id;
