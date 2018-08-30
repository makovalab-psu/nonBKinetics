set -e
set -x

while read folder; 
	do echo $folder; 
	mkdir $folder;
	cd $folder;

	#download all compressed files in a folder
	wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/hdf5/${folder}/*tar.gz
	wait

	for file in *tar.gz; do
		tar xvzf $file --xform='s#^.+/##x'; #extract compressed file
	done
	cd ..;
done < $1

