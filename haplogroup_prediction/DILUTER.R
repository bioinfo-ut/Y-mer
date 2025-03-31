# Simulating sequencing depth in k-mer values. 
# Input is counts file, format may differ

# Madala katvuse simuleerimine
# 

# Lubatud sisendid:

# Variant 1 (ilma päise ja k-meerita):

# 34
# 12
# ...


# Variant 2 (ilma päise aga k-meeriga):

# ATTTGAGGGGAGGGGAGG	34
# TTTTGAGGGGAGGGGACC	12
# ...


# Variant 3 (lühike päis koos k-meeriga):

#				Proov1	Proov2
# ATTTGAGGGGAGGGGAGG	34		23
# TTTTGAGGGGAGGGGACC	12		0
# ...


# Variant 4 (lühike päis ilma k-meerita):

#	Proov1	Proov2
#	34		23
#	12		0
# ...

# Variant 5 (Pikk päis k-meeriga või ilma - põhimõtteliselt mudeli loomiseks kasutatav fail):

#				Proov1	Proov2
#				A		A
# ATTTGAGGGGAGGGGAGG	34		23
# TTTTGAGGGGAGGGGACC	12		0
# ...


# Milliseid katvuseid tekitada (lisaks pannakse ka originaal faili kaasa, originaali katvusest suuremaid katvuseid ei soovita kasutada):
soovidevek=c(1, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001)


args = commandArgs(trailingOnly=TRUE)
algandmed = args[1]
valjundfail  = args[2]

# algandmed="C:/Maido/centro/NIPT/HalbJaHea/mudel/uus_nimekiri_NIPTIGA.txt"
# algandmed="C:/Maido/centro/NIPT/KatsedVer5/vanad.txt"

# valjundfail="C:/Maido/centro/NIPT/HalbJaHea/sim3.txt"


# Andmete sisselugemine
and0=read.table(algandmed, header=F, nrows=1, stringsAsFactors = FALSE)

if (prod(sapply(and0, is.character))==1) { 
	print("Sample ID-s on a first row")
	and00=read.table(algandmed, header=F, skip=1, nrows=1, stringsAsFactors = FALSE)
	if (prod(sapply(and00, is.character))==1) { 
		print("Ignoring 2. row") 
		and=read.table(algandmed, skip=2, header=FALSE) 
	    } else {
			and=read.table(algandmed, skip=1, header=FALSE) 
      }
   } else {
	print("No sample ID-s found")
	and=read.table(algandmed, header=F)
}

if (is.character(and[1,1])) {
   print("Ignoring 1. column - assuming it contains k-mers")
   and=and[,-1, drop=FALSE]
}

if (prod(sapply(and0, is.character))==1) {  colnames(and)=unlist(and0) }

katvused=colMeans(and)
katvused




# algandmed="C:/Maido/centro/NIPT/HalbJaHea/mudel/uus_nimekiri_NIPTIGA.txt"
# and=read.table(algandmed, skip=2, header=F, stringsAsFactors = FALSE)

# and=read.table("C:/Maido/centro/NIPT/HalbJaHea/mudel/HG00101_counts.txt", header=FALSE)

katseid_per_proov=(length(soovidevek)+1)
nridu=length(and[,1])
gc()
mat=matrix(0, nrow=nridu, ncol=katseid_per_proov*length(katvused))


for (j in 1:length(katvused)){	# Tsükkel üle proovide
	# Algne proov ise
	mat[,(j-1)*katseid_per_proov+1]=and[,j]

	for (i in 1:length(soovidevek)){ # Tsükkel üle soovitud katvuste
	  	mat[,(j-1)*katseid_per_proov+i+1] = rpois(nridu, soovidevek[i]/katvused[j]*and[,j])
	}
	
} # Tsükkel üle proovide


snimi=c("raw", paste(soovidevek, c("","A")[duplicated(soovidevek)+1], sep=""))
snimi=paste(snimi, c("","A")[duplicated(snimi)+1], sep="")

nimed=rep(snimi, length(katvused))
#nimed2=paste(colnames(and[j]), nimed, sep="_")
nimed2=paste(rep(colnames(and), each=length(snimi)), nimed, sep="_")


write.table(as.data.frame(t(nimed2)), file=valjundfail, col.names=FALSE, row.names=FALSE)
for (i in 1:(nridu/100000+99999/100000)){
  write.table(mat[((i-1)*100000+1):(min(i*100000, nridu)),], file=valjundfail, col.names=FALSE, row.names=FALSE, append=TRUE)
}


