#!/bin/bash
set -euo pipefail

#usage: ./get_one_ENA.sh <ftp_path> <outdir>

FTP=$1
OUT=$2

NAME=$(basename $FTP)

if [ ! -e "${OUT}/${NAME}" ]; then
	wget -O "${OUT}/${NAME}" ${FTP}
fi

