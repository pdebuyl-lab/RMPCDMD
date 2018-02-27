#!/bin/bash

filename="$1"

START=$(git log --follow --date=format:%Y --format="%ad" "${filename}" | tail -n 1)
END=$(git log --follow --date=format:%Y -n1 --format="%ad" "${filename}")

if [ "$START" = "$END" ] ; then
    YEARS="${START}"
else
    YEARS="${START}-${END}"
fi

# Keep only second line with '2c' to update existing headers
sed -i -e "1i ! This file is part of RMPCDMD" "${filename}"
sed -i -e "2i ! Copyright (c) $YEARS Pierre de Buyl and contributors" "${filename}"
sed -i -e "3i ! License: BSD 3-clause (see file LICENSE)\n"  "${filename}"

