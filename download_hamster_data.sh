#!/bin/bash

doy=265
hamster_dir="./C/lib/HAMSTER"

mkdir -p "$hamster_dir"

file=$(printf "DOY%03d.nc" "$doy")
url=$(curl -fsSL "https://opendata.physik.lmu.de/04zd8-7et52/?list" | grep "$file" | head -n 1)

wget -c "$url" -O "$hamster_dir/$file"