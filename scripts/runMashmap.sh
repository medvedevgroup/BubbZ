#Parameters description
#1) No. of genomes
#2) Path to the directory containing genomes (they should not contain anything else)

dir=out_$1
rm -rf $dir
mkdir $dir

rt=$dir/rt$1.txt

files=($2/*.fna)

cnt=0
end=$(($1 - 1))


declare -a input

for i in $(seq 0 $end);
do
	r=${files[$i]}
       	s="./mashmap -s 500 -r $r -q $r -t 1 -o $dir/mashmap_$i\_$i.tmp -f none"
	input+=("$s")
       	for j in $(seq $(($i + 1)) $end)
      	do
		q=${files[$j]}
		s="./mashmap -s 500 -r $r -q $q -t 1 -o $dir/mashmap_$i\_$j.tmp -f one-to-one"
		input+=("$s")
	done
done

#printf '%s\n' "${input[@]}"

align()
{
	bash -c "$1"
}

export -f align

printf '%s\n' "${input[@]}" | xargs -I @ -P 24 bash -c "align '@' "




