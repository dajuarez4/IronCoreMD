#!/bin/bash

MAIN_PATH="./"

echo "Checking QE jobs under:"
echo "$MAIN_PATH"
echo

DONE_COUNT=0
NOT_DONE_COUNT=0
TOTAL_COUNT=0

printf "%-12s  %-10s  %s\n" "STATUS" "FILE" "FOLDER"
printf "%-12s  %-10s  %s\n" "------" "----" "------"

while IFS= read -r outfile; do
    TOTAL_COUNT=$((TOTAL_COUNT + 1))
    folder=$(dirname "$outfile")
    file=$(basename "$outfile")

    if grep -q "JOB DONE." "$outfile"; then
        status="DONE"
        DONE_COUNT=$((DONE_COUNT + 1))
    else
        status="NOT_DONE"
        NOT_DONE_COUNT=$((NOT_DONE_COUNT + 1))
    fi

    printf "%-12s  %-10s  %s\n" "$status" "$file" "$folder"

done < <(find "$MAIN_PATH" -type f -name "*fe1.out" | sort)

echo
echo "Summary:"
echo "Total fe1.out files found : $TOTAL_COUNT"
echo "DONE                     : $DONE_COUNT"
echo "NOT DONE                 : $NOT_DONE_COUNT"
