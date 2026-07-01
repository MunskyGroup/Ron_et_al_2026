#!/usr/bin/env bash

file="DUSP1_SSITcellresults_Final_Sep18.csv"
colname="replica"

awk -F',' -v col="$colname" '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        gsub(/^"|"$/, "", $i)
        if ($i == col) {
            c = i
            break
        }
    }
    if (!c) {
        print "Column \"" col "\" not found." > "/dev/stderr"
        exit 1
    }
    next
}

NR > 1 {
    val = $c
    gsub(/^"|"$/, "", val)

    if (val != "") {
        if (!seen[val]++) {
            values[++n] = val
        }
    }
}

END {
    if (n == 0) {
        print "No values found in column \"" col "\"." > "/dev/stderr"
        exit 1
    }

    print "Values in column \"" col "\":"
    for (i = 1; i <= n; i++) {
        print values[i]
    }
}
' "$file"
