#!/usr/bin/bash

cat filename | while read line

do
find-peaks ./chunks/$line cons.schema AGO-Signal
done
