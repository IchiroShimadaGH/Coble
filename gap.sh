#!/bin/sh



if [ $# -eq 0 ] ; then
  mem=4G
else
  mem=$1
fi;

echo "memory=" "$mem"

/usr/local/Cellar/gap/4.15.1/bin/gap -L ./gapReads/startgwsp/startws.ws -o $mem
