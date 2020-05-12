#Parameters description
#1) No. of genomes
#2) Path to the directory containing genomes (they should not contain anything else)

echoerr() { echo "$@" 1>&2; }

dir=out_$1
rm -rf $dir
mkdir $dir

rt=$dir/rt$1.txt

files=($2/*.fna)
cnt=0
end=$(($1 - 1))
for i in $(seq 0 $end);
do
        r=${files[$i]}
        cat $r >> out_$1/genomes.fa
done

k=21
d=$1
outdir=out_$d
gfile=out_$1/genomes.fa
dbfile=$outdir/salmonella_$d\_$k.dbg
/usr/bin/time -f "TwoPaco: %e elapsed %M memory" ./twopaco -t 16 -k $k -f 37 $gfile -o $dbfile -a 9600 2> $outdir/rt.txt
/usr/bin/time -f "LCB: %e elapsed %M memory" ./bubbz-map $gfile --graph $dbfile -k $k -b 200 -o $outdir -m 250 -t 24 > $outdir/log.txt 2>> $outdir/rt.txt




