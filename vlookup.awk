FNR==NR {
    a[$1]=$2
    next
}
{ if ($2 in a) {printf "%s\t%s\t%s\n",$1,$2,a[$2]} else {printf "%s\t%s\tNA\n",$1,$2} }
