#!/bin/bash
# DC script to set up GVA analyses, one script per gene

geneList=/home/rejudcu/reference/allGenes101115.txt
# geneList=/home/rejudcu/reference/DRDgenes.txt

# disease="UCLEx.Prionb2"
# model="ExAC.ct08.rare"
# model="ct08.cleaned"
# must be in this order or else qdel will delete all the ct08 jobs
disease="SSS2"
model="vrare.1"

homeFolder=/cluster/project9/bipolargenomes/scores
argFolder=/home/rejudcu/pars
softwareFolder=/home/rejudcu/bin
dataHome=/home/rejudcu
copyVCF=yes
# only attempt to copy vcf files if not too big, else arg file must specify absolute path to vcf files

if [ -z "$disease" -o -z "$model" ]
then
	echo Error in $0: must set environment variables disease and model
	exit
fi

someGenesLeft=no

for d in $disease
do
for m in $model
do

testName=$d.$m
argFile=$argFolder/gvs.$testName.arg


# workFolder=/cluster/project8/bipolargenomes/GVA

workFolder=$homeFolder/$d/$testName
mkdir $homeFolder/$d
mkdir $workFolder

nSplits=100

splitScript=$workFolder/scripts/split${nSplits}s.sh
scriptName=$testName.runSplit${nSplits}.sh
mainSplitScript=$workFolder/scripts/$scriptName

qdel $testName.'runSplit*'

nhours=4
vmem=6 
memory=2
queue=queue6
scratch=10

if [ ! -e $workFolder ]; then mkdir $workFolder; fi;
wastebin=$workFolder/wastebin
if [ ! -e $wastebin ]; then mkdir $wastebin; fi
if [ ! -e $workFolder/results ]; then mkdir $workFolder/results; fi;
if [ -e $workFolder/error ]; then mv $workFolder/error $wastebin/error; ( rm -r $wastebin/error & ) ; fi;
mkdir $workFolder/error
if [ -e $workFolder/scripts ]; then mv $workFolder/scripts $wastebin/scripts; (rm -r $wastebin/scripts & ); fi;
mkdir $workFolder/scripts; 
if [ -e $workFolder/temp ]; then mv $workFolder/temp $wastebin/temp; (rm -r $wastebin/temp & ); fi;
mkdir $workFolder/temp; 

cat $geneList | while read geneName
    do
    shellScript=$workFolder/scripts/runGVS.$testName.$geneName.sh
    if [ -e $shellScript ] ; then rm $shellScript; fi
    outFile=$workFolder/results/$testName.$geneName.sao
	scoreFile=$workFolder/results/$testName.$geneName.vsco
	elogFile=$workFolder/results/$testName.$geneName.elog
# I am going to add an exclusion log file so I can find which variants failed which conditions
    if [ ! -e $outFile ]
    then 
		echo "PATH=$softwareFolder:\$PATH 
		tempFolder=$geneName
		mkdir \$tempFolder 
		cd \$tempFolder 
		rm -f gva*.$geneName.* gvs*.$geneName.*
		commLine=\"geneVarAssoc --arg-file $argFile --gene $geneName \" 
		echo Running:
		echo \$commLine
		\$commLine 
		echo finished running geneVarAssoc
		if [ -e gva.$geneName.dat ]
		then
			commLine=\"getVarScores --flagfile ../reference/varFlags.txt --gcdatafile gva.$geneName.dat --locusfilterfile gva.$geneName.lf.par --locusweightfile gva.$geneName.lw.par --locusnamefile gva.$geneName.comm.par --filterfile gva.$geneName.filter.par --outfile gvs.$geneName.sao --scorefile gvs.$geneName.vsco --weightfactor 1.000000 \"
			echo Running:
			echo \$commLine
			\$commLine 
			echo finished running getVarScores
		else
			echo geneVarAssoc did not extract valid variants for $geneName
			echo geneVarAssoc did not extract valid variants for $geneName > $outFile
		fi
		# cp gvs.$testName.$geneName.elog $elogFile
		cp gvs.$geneName.vsco $scoreFile 
		cp gvs.$geneName.sao $outFile 
		cd ..
		# avoid upsetting bash by removing the current working directory
		rm -r \$tempFolder" >> $shellScript
    fi
    done

# was \$commLine > gvs.$testName.$geneName.elog 
		
nScriptsWritten=`find $workFolder/scripts -name 'runGVS.*.sh' | wc -l`
if [ $nScriptsWritten -lt $nSplits ]
then 
	nSplits=$nScriptsWritten
fi 

if [ -e  $mainSplitScript ] ; then rm  $mainSplitScript; fi

echo "
#!/bin/bash
#$ -S /bin/bash
#$ -e $workFolder/error
#$ -o $workFolder/error
#$ -l tscr=${scratch}G
#$ -l tmem=${vmem}G,h_vmem=${vmem}G
#$ -l h_rt=${nhours}:0:0
#$ -t 1-$nSplits
#$ -V
#$ -R y

# I used to have #$ -cwd but I am going to try just omitting it as sometimes cannot cd to it
# If that does not work may try -wd /scratch0

date
echo bash -x $splitScript \$SGE_TASK_ID
bash -x $splitScript \$SGE_TASK_ID
date
" > $mainSplitScript

echo "
#!/bin/bash
set +e
#  was exiting after running just one, possibly because no proper exit code from script
# this should switch off errexit
echo Running \$0 with argument \$1
mkdir /scratch0/$USER
tmpDir=/scratch0/$USER/\$RANDOM
mkdir \$tmpDir
cd \$tmpDir
if [ $copyVCF = yes ]
then
	mkdir vcf
	mkdir vcf/$disease
	cp -L $dataHome/vcf/$disease/* vcf/$disease
	# copy destination of links
fi
mkdir reference
cp -L $dataHome/reference/* reference
mkdir pars
cp -L $dataHome/pars/* pars
# runGVA.sh will create a temporary folder named after the gene and cd to it
n=1
find $workFolder/scripts -name 'runGVS*sh' | while read f
do
if [ .\$n == .\$1 ]
then
	echo running source \$f # try using source $f instead of bash $f
	source \$f
	echo finished running source \$f
fi
if [ \$n -eq $nSplits ]
then
	n=1
else
	n=\$(( \$n + 1 ))
fi
done
rm -r \$tmpDir
" > $splitScript

count=`find $workFolder/scripts -name 'runGVS*sh' | wc -l`

if [ $count -gt 0 ]
then
	echo wrote $count scripts
	echo qsub -N $scriptName $mainSplitScript
	pushd $workFolder
# reason for this is that I would get Eqw with qstat -j error message: error: can't chdir to /home/rejudcu/tmp: No such file or directory 
	qsub -N $scriptName $mainSplitScript
	popd
	someGenesLeft=yes
else
	echo No genes left to do for $testName	
fi

done
done

logFile=${0##*/}
if [ $someGenesLeft = yes ]
then
	echo will schedule script to run again
	echo "export disease=\"$disease\"; export model=\"$model\"; bash $0 &> $workFolder/$logFile.log" | at now + $nhours hours
else
	echo date > $workFolder/$logFile.log
	echo All results files written OK >> $workFolder/$logFile.log
fi
