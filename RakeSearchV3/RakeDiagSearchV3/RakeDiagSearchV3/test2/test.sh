#!/bin/bash

rm -f boinc_finish_called checkpoint.txt result.txt stderr.txt

time ./rakesearch10 > /dev/null
#time ../tests/sde-external-8.12.0-2017-10-23-win/sde -skx -- ./rakesearch10

diff -sq result.txt result.txt.ref
