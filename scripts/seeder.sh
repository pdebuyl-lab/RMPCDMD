#!/bin/bash

SEED=$(dd count=1 bs=8 if=/dev/urandom 2>/dev/null | od -A n -t d8)

if [ ! -n "$1" ]
then
echo ${SEED}
exit
fi  

sed -e "s/\${SEED}/${SEED}/" "$1"
