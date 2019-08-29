# Annotation pre-processing step, modified to accommodate for the ">" in front of my fasta lines.
fasta="STAR_fasta_tr.fa"
gtf="Homo_sapiens.GRCh38.92.chr.gtf"

outputFile <- "out_ex"

system(paste('grep ">" ',fasta,' > fastanames.txt',sep=''))#'grep ">" /mnt/mr01-data01/mqbssppe/cufflinksHumanAnnotation/transcriptome_data/known.fa > fastanames.txt')
myFile <- file("fastanames.txt",open="r")
K <- as.numeric(strsplit(system("wc -l fastanames.txt",intern = TRUE),split = " ")[[1]][1])
trNames <- numeric(K)
i <- 0
while (length(oneLine <- readLines(myFile, n = 1, warn = FALSE)) > 0){
  i <- i + 1
  x = strsplit(oneLine,split = " ")[[1]][1]
  trNames[i] <- substring(x, 2, nchar(x))
}
mean(is.na(trNames))
head(trNames); tail(trNames)
close(myFile)
write.table(trNames,file = "transcriptNames.tr",quote = FALSE,row.names = FALSE, col.names = FALSE)
#################################################################################################################

myFile <- file(gtf,open = "r")
d1 <- 100000
temp <- array(data = NA,dim = c(d1,2))
i <- 0
while(length(line <- readLines(myFile, n = 1, warn = FALSE)) > 0){
  # print(line)
  i <- i + 1
  ln1 <- strsplit(line,split=c("gene_id"))[[1]][2]
  ln2 <- strsplit(ln1,split = "transcript_id")
  geneName <- strsplit(ln2[[1]][1],split= "\"")[[1]][2]
  trName <- strsplit(ln2[[1]][2],split = "\"")[[1]][2]
  if(i > d1){
    cat(paste("line:",i),"\n")
    d1 <- d1 + 100000
    tempNew <- array(data = NA,dim = c(d1,2))
    tempNew[1:(i-1),] <- temp
    temp <- tempNew
    tempNew <- 0
  }
  temp[i,] <- c(geneName,trName)
}
# Some NAs at the beginning...
temp <- temp[1:i,]
mean(is.na(temp))
head(temp); tail(temp)

close(myFile)

myNames <- unique(temp)
if(length(trNames)!= dim(myNames)[1]){message("number of transcripts in fasta file is not equal to number of transcripts in gtf file")}else{message("number of transcripts in fasta file is equal to number of transcripts in gtf file")}
notFoundTranscripts <- which(is.na(match(trNames,myNames[,2]))==TRUE)
if (length(notFoundTranscripts) > 0){message("not compatible transcripts detected")}
colnames(myNames) <- c("geneID","transcriptID")

head(myNames); tail(myNames)

write.table(myNames,file = "geneANDtranscripts.txt",quote =FALSE)

####################################################################################################################
l <- K
finalNames <- trNames

knownIsoformNames <- myNames
v <- knownIsoformNames
knownIsoformNames <- cbind(as.character(v[,1]),as.character(v[,2]))
noNameIndex <- which(knownIsoformNames[,2] == "")
knownIsoformNames[noNameIndex,2] <- knownIsoformNames[noNameIndex[1],1]

perm <- match(finalNames,as.character(knownIsoformNames[,2]))

knownIsoformNames <- knownIsoformNames[perm,]
knownIsoformNames <- cbind(knownIsoformNames,1:l)

nTr <- numeric(l)
for(i in 1:l){
  j = knownIsoformNames[i,1]
  nTr[i] <- length(which(knownIsoformNames[,1]==j))
}
knownIsoformNames <- cbind(knownIsoformNames,nTr)
trIDperGene <- numeric(l)
tt <- table(knownIsoformNames[,1])
gMax <- length(tt)
for(i in names(tt)){
  myIndex <- which(knownIsoformNames[,1] == i)
  trIDperGene[myIndex] <- 1:length(myIndex)
}

knownIsoformNames <- cbind(knownIsoformNames,trIDperGene)

dim(knownIsoformNames)
head(knownIsoformNames); tail(knownIsoformNames)

write.table(knownIsoformNames,outputFile,row.names=FALSE,col.names=c("geneID","trName","trID","nTrPerGene","trIDperGene"),quote=FALSE)
