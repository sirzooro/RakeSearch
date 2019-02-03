#!/bin/sh

rm -f boinc_finish_called checkpoint.txt result.txt
echo Started RakeSearch test...
time ./rakesearch
diff -qs result.txt result-ok.txt

