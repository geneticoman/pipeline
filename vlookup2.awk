FNR==NR {
    a[$1]=$3
    next
}
{ if ($1 in a) {printf "%s\t%f\n",$1,a[$1]} else {printf "%s\t0.000000\n",$1} }
