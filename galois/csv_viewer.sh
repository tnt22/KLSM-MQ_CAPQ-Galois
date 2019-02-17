# Usage csv_viewer.sh results.csv
# 
# This tool allow you to easily view csv files in the console. Allowing you to make a comparation between the columns and rows.
column -s, -t < $1 | less -#2 -N -S

