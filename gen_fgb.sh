#rm anfr.fgb
#rm anfr.json anfr2.json
#FlatGeobuf
# ogr2ogr -f GeoJSON anfr_all.json anfr2026.gpkg -sql "select sites.id, CAST(sites.id AS INTEGER)+0 as site_id, sites.dms, lon, lat, address, postcode, inseecode, geom,
#     count(distinct operator) operator_count, cast(group_concat(distinct operator) as text) operator_list,
#     count(distinct supports.id) support_count, max(sup_height) h_max,
#     count(distinct station_id) sta_count,
#     count(distinct antennas.id) ant_count,
#     count(distinct system) tech_count, cast(group_concat(distinct system) as text) tech_list,
#     count(distinct transmitters.id) tx_count,
#     count(distinct bandstr) band_count,
#     sites.tech_bitmask1, sites.tech_bitmask2
#     from sites
#     inner join antennas on antennas.sup_id=supports.id
#     inner join supports on sites.id=supports.site_id
#     inner join transmitters on transmitters.antenna_id=antennas.id
#     inner join id_systems on transmitters.system_id=id_systems.id
#     inner join v_bandgroups on transmitters.bandgroup_id=v_bandgroups.id
#     inner join stations on stations.id=transmitters.station_id
#     inner join id_operators on id_operators.id=stations.operator_id
#     group by sites.id"

#ogr2ogr -f GeoJSON anfr2.json anfr2026.gpkg -sql "select sites.id site_id, geom from sites"

rm anfr.json anfr.json.gz
ogr2ogr -f GeoJSON anfr.json anfr2026.gpkg -sql "select sites.id, CAST(sites.id AS INTEGER)+0 as site_id, geom,
    count(distinct supports.id) support_count, max(sup_height) h,
    count(distinct operator) operator_count,
    cast(group_concat(distinct operator) as text) operator_list,
    count(distinct station_id) sta_count,count(distinct antennas.id) ant_count,
    sites.tech_bitmask1, sites.tech_bitmask2
    from sites
    inner join antennas on antennas.sup_id=supports.id
    inner join supports on sites.id=supports.site_id
    inner join transmitters on transmitters.antenna_id=antennas.id
    inner join id_systems on transmitters.system_id=id_systems.id
    inner join stations on stations.id=transmitters.station_id
    inner join id_operators on id_operators.id=stations.operator_id
    group by sites.id"
gzip anfr.json

sqlite3 -json anfr2026.gpkg "select min(support_count) sup_c_min, max(support_count) sup_c_max, min(h) h_min, max(h) h_max,
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
    group by sites.id)"


