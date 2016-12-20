#Let the field seperator to a newline
IFS="
"
# Loop through the file
for line in `cat $1`;do
# Echo the line (echo could be changed to whatever command you want)
    grep -v $line $2 > $3
    cp $3 $2
done
