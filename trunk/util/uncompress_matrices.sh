#! /bin/sh

find Matrices -name '*.sms.gz' -print0 \
    | xargs --null -n1 gunzip -v
