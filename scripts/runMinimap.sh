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
       	s="./minimap2 -DP -k19 -w19 -m200 -t1 $r $r > $dir/minimap_$i\_$i.tmp"
	input+=("$s")
       	for j in $(seq $(($i + 1)) $end)
      	do
		q=${files[$j]}
		s="./minimap2 -D -cx asm10 -t1 $r $q > $dir/minimap_$i\_$j.tmp"
		input+=("$s")
	done
done

align()
{
	/usr/bin/time -f "Minimap: %e elapsed %M memory" bash -c "$1"
}

export -f align

printf '%s\n' "${input[@]}" | xargs -I @ -P 24 bash -c "align '@' "



