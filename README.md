# fradiodb
This tool helps process the [French open data on radio sites > 5W](https://www.data.gouv.fr/fr/datasets/donnees-sur-les-installations-radioelectriques-de-plus-de-5-watts-1/).
That data actually has a some duplications, and the tables and field names as well as the schema are not always very straightfoward, therefore some processing is done in order to obtain the following target schema:
![schema](schema.png)

It is to be noted that
* Table `sites` was created because there may be several supports at the same location
* One support obviously has several stations, but conversely the current definition of "station" may encompass several transmitters that (due to historical reasons) may be on different supports
* `bandgroups.id` is not a PK (it's not unique). I could have defined another table "bands_id" where bands would have been unique, but this would add complexity and overall the amount of band combinations is only 2x the amount of distinct bands

Then, the tool also enable the creation of tables `gen_sites` and `gen_sectors` (possibly in another DB) to summarize within a single table some useful information from the standpoint of sites (i.e. geographical address) and sectors (i.e. with a boresight), e.g. for the purpose of importing in a GeoDataFrame and mapping with leaflet.js.

## Warning: work-in-progress
This is a work-in-progress, and both the schema and the processing may still have bugs and/or change
