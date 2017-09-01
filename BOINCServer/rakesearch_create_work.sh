#!/bin/sh

# Workunit creation for the specified rank and number of a square

IN_FILE="/home/boincadm/sandbox/files/wu_000014.txt"
RANK="8"

BOINC_INFILE=`basename $IN_FILE`
IN_NAME=`basename $IN_FILE .txt`
WU_N=`echo $IN_NAME | cut -d'_' -f2` 
WU_NAME=${RANK}_${WU_N}

./bin/stage_file --copy ${IN_FILE}
./bin/create_work --appname rakesearch --wu_name ${WU_NAME} ${BOINC_INFILE}

