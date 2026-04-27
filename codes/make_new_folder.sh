#!/bin/bash

OLD_TEMP="5000K"
NEW_TEMP="4500K"

OLD_NUM="5000"
NEW_NUM="4500"

mkdir -p "$NEW_TEMP"

for a in $(seq -f "%.2f" 2.29 0.01 2.49); do
    src="${a}_${OLD_TEMP}"
    dst="${NEW_TEMP}/${a}_${NEW_TEMP}"

    if [[ ! -d "$src" ]]; then
        echo "[SKIP] Missing source folder: $src"
        continue
    fi

    mkdir -p "$dst"

    echo "[COPY] $src -> $dst"

    # Copy files
    cp -v "$src"/fe*.in "$dst"/ 2>/dev/null
    cp -v "$src"/run_jakar.sbatch "$dst"/ 2>/dev/null
    cp -v "$src"/live* "$dst"/ 2>/dev/null

    # Change folder references like 2.29_5000K -> 2.29_4500K
    sed -i "s/${a}_${OLD_TEMP}/${a}_${NEW_TEMP}/g" "$dst"/* 2>/dev/null

    # Change 5000K -> 4500K inside copied files
    sed -i "s/${OLD_TEMP}/${NEW_TEMP}/g" "$dst"/* 2>/dev/null

    # Change 5000 -> 4500 inside fe input and sbatch file
    sed -i "s/${OLD_NUM}/${NEW_NUM}/g" "$dst"/fe*.in 2>/dev/null
    sed -i "s/${OLD_NUM}/${NEW_NUM}/g" "$dst"/run_jakar.sbatch 2>/dev/null

    # Change pseudo_dir path because files are now one folder deeper
    sed -i "s@pseudo_dir *= *'../../pseudo'@pseudo_dir = '../../../pseudo'@g" "$dst"/fe*.in 2>/dev/null

done

echo "Done creating folders inside $NEW_TEMP/"