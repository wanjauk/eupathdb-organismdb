#!/usr/bin/env bash
#
# TriTrypDB organism package generation post-processing
# Keith Hughitt (khughitt@umd.edu)
# 2015/01/29
#
# This package applies some post-processing to the output of the previous R
# scripts, renaming packages and fixing metadata fields.
##############################################################################

# TODO: Add TriTrypDB DBSCHEMA?
# https://github.com/genome-vendor/r-bioc-annotationdbi/tree/master/inst/DBschemas

# YAML parser
# http://stackoverflow.com/questions/5014632/how-can-i-parse-a-yaml-file-from-a-linux-shell-script 
function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

# load configuration
eval $(parse_yaml config.yaml)

# old package name (e.g. "org.Lmajor.eg.db")
printf -v orgdb_name_old 'org.%s%s.eg.db' ${genus:0:1} ${species}

# name without the .db suffix
orgdb_name_short_old=${orgdb_name_old/.db/}
orgdb_name_short=${orgdb_name/.db/}

# rename and enter directory
mv $orgdb_name_old $orgdb_name
cd $orgdb_name

# Fix DESCRIPTION
sed -i "s/$orgdb_name_old/$orgdb_name/" DESCRIPTION
sed -i "s/species:.*/species: $description/" DESCRIPTION
sed -i "s/Entrez/TriTrypDB/" DESCRIPTION

# Fix NAMESPACE
sed -i "s/$orgdb_name_short_old/$orgdb_name_short/" NAMESPACE

# Fix sqlite database
dbpath=inst/extdata/${orgdb_name_short}.sqlite
mv inst/extdata/${orgdb_name_short_old}.sqlite $dbpath

chmod +w $dbpath
sqlite3 $dbpath "UPDATE metadata SET value=\"$description\" WHERE name='SPECIES';"
chmod -w $dbpath

# Fix manual pages
for suffix in "BASE.Rd" "ORGANISM.Rd" "_dbconn.Rd"; do
    mv man/${orgdb_name_short_old}${suffix} man/${orgdb_name_short}${suffix}
    sed -i "s/$orgdb_name_short_old/$orgdb_name_short/" man/${orgdb_name_short}${suffix} 
done

# Fix zzz.R
sed -i "s/$orgdb_name_short_old/$orgdb_name_short/" R/zzz.R

echo "Done!"
