#!/bin/bash

export DDP=1
export FIELD='G9'
export GOLD_DIR=~/data/GAMA4/

if [ -z "$DDP" ]; then
    echo 'Skipping DDP randoms.'
else
    echo 'Generating DDP randoms.'

    DDP_FILE=$GOLD_DIR/gama_gold_ddp.fits

    if [ -f "$DDP_FILE" ]; then
        echo "Found DDP file"
    else       
	echo $DDP_FILE' does not exist.'
        exit 1
    fi
    
    ddp1_zmin=$(fitsheader $DDP_FILE -k 'DDP1_ZMIN' -e 1)
    ddp1_zmin=$(echo $ddp1_zmin | cut -d " " -f 9)

    ddp1_zmax=$(fitsheader $DDP_FILE -k 'DDP1_ZMAX' -e 1)
    ddp1_zmax=$(echo $ddp1_zmax | cut -d " " -f 9)

    echo $ddp1_zmin
    echo $ddp1_zmax

    python ~/DESI/randoms.py --field $FIELD --prefix randoms_ddp1 --zmin $ddp1_zmin --zmax $ddp1_zmax $DRYRUN > $CODE_ROOT/logs/randoms_ddp_$FIELD.log
fi
