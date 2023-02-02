#!/bin/bash

cutoff_energy=-4.0

rm box*.txt AllSets.pdb AllScores.txt config.txt PredictedInternal.pdb temp_PredictedWaters*.pdb

ending="qt"
pdbqt=$1$ending

if [[ ! -s $pdbqt ]]; then
 ./pythonsh prepare_receptor4.py -r $1 -A bonds_hydrogens
fi
if [[ ! -s $pdbqt ]]; then
 ./pythonsh prepare_receptor4.py -r $1 -A hydrogens
fi
if [[ ! -s $pdbqt ]]; then
 ./pythonsh prepare_receptor4.py -r $1
fi

./pdbqt $pdbqt
cp new.pdbqt $pdbqt
rm new.pdbqt

#./make_config $pdbqt

#box_size=$2

#if [[ $2<0 ]]; then
# box_size="10"
#fi

#./conf $box_size

# read box number into the variable num_boxes
#read num_boxes<box_number.txt
num_boxes=$(ls ../boxes/box* -1q | wc -l)
echo "There are $num_boxes tunnel points"
max_cluster=3

for kconf in $(seq 1 $num_boxes); do
echo working on box_$kconf.txt...
for j in $(seq 1 $max_cluster); do
	progress1="Processing part "$kconf" out of "$num_boxes
	progress2="Clustering round "$j" out of "$max_cluster
	echo ""
	echo "Dowser with Caver, software for hydrating tunnels in protein structures"
	echo ""
	echo "By Panyue Wang at University of California, Davis"
	echo ""
	echo ""
	echo $progress1
	echo ""
	echo $progress2
	echo ""
	./run --receptor $pdbqt --ligand water.pdbqt --config ../boxes/box_"$kconf".txt  --out Set"$j".pdbqt --exhaustiveness=20 > /dev/null
	grep OW Set"$j".pdbqt >> AllSets.pdb
	grep RESULT Set"$j".pdbqt | awk '{print $4}' >> AllScores.txt 
	rm Set"$j".pdbqt
done

echo load_data
./load_data
rm AllScores.txt AllSets.pdb

R --slave --vanilla --quiet --no-save  <<RSCRIPT

#### Loading and preparing coordinate data.
	pdb <- read.table("data.pdb")		# Loading the table of docking results
	coords <- pdb[,1:4]				# Using the coordinates of the PDB file or...

#### First round of clustering. 
	d <- dist(coords)				# The matrix of euclidean distances
	fit <- hclust(d, method="single") 		# Ward minimum variance clustering.
	clusts <- cutree(fit,h=0.5)			# h is the distance cut off. 

	#plot(fit)					# In case you want to visualise the tree...
	#rect.hclust(fit, h=1.6, border="red") 		# ... and view the clusters at the cutoff.

# Taking the mean position of the members of each cluster as the represenative.	
	clustcoords <- c()					
	for (j in 1:max(clusts)){
		question <- t(as.matrix(coords[clusts==j,])) 		 
		if (nrow(question) == 1){
			clustcoords <- rbind(clustcoords,round(question,digits=3))
		} else {
			clustcoords <- rbind(clustcoords,round(apply(coords[clusts==j,],2,mean),digits=3))
		}
	}

#### Second round of clustering.
	fit2 <- hclust(dist(clustcoords),method="single")
	clusts2 <- cutree(fit2,h=1.6)
	clustcoords2 <- c()					
	for (j in 1:max(clusts2)){
		question <- t(as.matrix(clustcoords[clusts2==j,])) 		
		if (nrow(question) == 1){
			clustcoords2 <- rbind(clustcoords2,round(question,digits=3))
		} else {
			clustcoords2 <- rbind(clustcoords2,round(apply(clustcoords[clusts2==j,],2,mean),digits=3))
		}
	}
	Coords <- as.matrix(clustcoords2)		# These are my water final coordinates.


#### Outputting the coordinates into the correct pdb format.
# Setting up the column distances.
	atom.print <- function(card = "ATOM", eleno, elety, alt = "",resid, chain = "", resno, insert = "", x, y, z, o = "1.00", b, segid = "") {
        	format <- "%-6s%5s  %-3s%1s%-4s%1s%4s%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s"
       		if (nchar(elety) > 3) {
			format <- "%-6s%5s %-4s%1s%-4s%1s%4s%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s"
        	}
        	sprintf(format, card, eleno, elety, alt, resid, chain, resno, insert, "", x, y, z, o, b, "", segid)
	}

# Pre-assignment for below.
	eleno <- 1:nrow(Coords)
	resno <- 1:nrow(Coords)
	chain <- rep("A",nrow(Coords))
	elety <- rep("O",nrow(Coords))
	resid <- rep("HOH",nrow(Coords))
	o <- rep(1.00,nrow(Coords))

# Putting together the coordinates in the correct format.
	tbl<-c()
	for (j in 1:nrow(Coords)) {
		coords <- matrix(round(as.numeric(Coords[j,]), 4), ncol = 4, byrow = TRUE)
	tbl <- rbind(tbl, atom.print(eleno = eleno[j], elety = elety[j], resid = resid[j], chain = chain[j], resno = resno[j], x = coords[1, 1], y = 		coords[1,2], z = coords[1, 3], o = o[j], b = coords[1, 4]))
	}

# Checking to see if the file name "temp_PredictedWaters.pdb" is already taken. If not, a number is added at the end.
	folderlist <- list.files()
	file <- "temp_PredictedWaters.pdb"	
	k <- 2
	while(is.element(file,folderlist)){
		file <- paste("temp_PredictedWaters_",k,".pdb",sep="")
		k <- k+1
	}

# Writing the header, pdb data and footer details into the file.
#	line1 <- "TITLE     Binding Site Water Molecules Predicted Using the WaterDock Method $kconf"
#	line2 <- paste("REMARK   ",date())
#	line3 <- "REMARK    Written by Gregory Ross"

#	write(line1,file=file,ncolumns=1000,append=FALSE)
#	write(line2,file=file,ncolumns=1000,append=TRUE)
#	write(line3,file=file,ncolumns=1000,append=TRUE)

	write.table(tbl, file = file, quote = FALSE, row.names = FALSE,col.names = FALSE, append = TRUE)
#	cat("TER   ", "END", file = file, sep = "\n", append = TRUE)

RSCRIPT

	rm data.pdb

done

#rm box*.txt # $pdbqt
cat temp_PredictedWaters*.pdb > PredictedWaters.pdb
rm temp_PredictedWaters*.pdb 

./reform -pdbin $1 -pdbout reform.pdb
cat reform.pdb PredictedWaters.pdb > reform2.pdb
#./drain reform2.pdb PredictedWaters.pdb surface.wat PredictedInternal.pdb
#rm surface.wat PredictedWaters.pdb 
#
#./sorting
#rm PredictedInternal.pdb
#./internal_predicted
#rm sorted.pdb

#./placeWat reform.pdb PredictedInternal.pdb rotate > placed.pdb 
#
#cp placed.pdb refined.pdb
#./choosing $cutoff_energy
#rm placed.pdb refined.pdb reform.pdb reform2.pdb



