#!/usr/bin/awk -f

# Example:
# translator.awk dict.txt input.txt > translated_input.txt

NR == FNR {
  rep[$1] = $2
  next
} 

{
  for (key in rep)
    gsub(key, rep[key])
  print
}