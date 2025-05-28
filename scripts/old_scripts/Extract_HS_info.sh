#!/bin/bash
SERVER=/run/user/1000/gvfs/smb-share:server=maximus-nas.local,share=genomics
BAM=/media/Data/WES/BAM/

find ${BAM}/*/ -type f -name "*_hs_metrics.txt" | while read file; do awk 'NR==8 {print FILENAME "\t" $0}' "$file"; done | awk -F '/' '{print $(NF-1), $0}' OFS='\t' | cut -f 1,3- > HS_info.txt

head $(ls -d /media/Data/WES/BAM/*/*_hs_metrics.txt | head -n 1) -n 7 | tail -n 1 > heading.txt

cat heading.txt HS_info.txt | awk 'NR==1 {print "Sample.ID\t" $0; next} 1'> ${BAM}/HS_combined.txt

rm HS_info.txt
rm heading.txt

