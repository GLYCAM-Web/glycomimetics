rm -f summary_csv.txt
cd all_csvs
for i in *.csv;do
    tail -n 1 ${i} | awk -v moiety="${i}" '{print moiety" "$2"\t"$3}' 
    tail -n 1 ${i} | awk -v moiety="${i}" '{print moiety" "$2"\t"$3}' >> summary_csv.txt
done
cd ..
