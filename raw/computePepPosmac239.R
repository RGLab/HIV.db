
#########################################################
##1.pair wise alignment of consecutive peptide sequence to construct entire clade sequence
#######################################################
res<-vector("list",2)
names(res)<-c("mac239", "E660")

#get peptides for E660
JPT239<-read.table("~/Programs/git/HIV.db/inst/extdata/JPTmac239.gal", header=TRUE)
JPTO<-JPT239[with(JPT239, order(Name)),]
strs<-strsplit(as.character(JPTO[,1]), split="-")
strs2<-unlist(lapply(strs, function(x){
			tail(x,1)
		}))

pepE660<-as.character(JPTO$ID[which(strs2!="239")])
pep239<-as.character(JPTO$ID[which(strs2!="E660")])

#E660
start<-end<-rep(0,length(pepE660)) #initialize start/end
tmpSeq<-do.call(function(x,...){paste(x,sep="",...)},as.list(rep("-",886)))
start[1]<-1
end[1]<-start[1]+nchar(pepE660[1])-1
substr(tmpSeq,start[1],end[1])<-pepE660[1]
for(i in 2:length(pepE660))
{
	tmp<-pairwiseAlignment(pepE660[i],pepE660[i-1],type="overlap")
	start[i]<-start[i-1]+start(tmp@subject@range)-1
	end[i]<-start[i]+nchar(pepE660[i])-1
	substr(tmpSeq,start[i],end[i])<-pepE660[i]
}
res[["E660"]]<-list(peptide=pepE660,start=start,end=end,consensus=substr(tmpSeq,start[1],end[length(pepE660)]))

#239
start<-end<-rep(0,length(pep239)) #initialize start/end
tmpSeq<-do.call(function(x,...){paste(x,sep="",...)},as.list(rep("-",879)))
start[1]<-1
end[1]<-start[1]+nchar(pep239[1])-1
substr(tmpSeq,start[1],end[1])<-pep239[1]
for(i in 2:length(pep239))
{
	tmp<-pairwiseAlignment(pep239[i],pep239[i-1],type="overlap")
	start[i]<-start[i-1]+start(tmp@subject@range)-1
	end[i]<-start[i]+nchar(pep239[i])-1
	substr(tmpSeq,start[i],end[i])<-pep239[i]
}
res[["mac239"]]<-list(peptide=pep239,start=start,end=end,consensus=substr(tmpSeq,start[1],end[length(pep239)]))

SVres<-res

#########################################################
#2.using MUSCLE to perform alignment of multipe sequence
#....................
#######################################################

######################################################################
##3.compute the new HXB2 positions for peptides based on multi-aligned clade sequence
###################################################################


######################################################################
##extended substr function to extract a cernain length of substring 
#without couting ignored sybmols
######################################################################

substrEx<-function(x,start,len,ignore="-"){
	res<-NULL
	letterCount<-0##counter of the number of amino acides
	while(start<=nchar(x))##sequentially read letter by letter within current aligned clade sequence
	{
		curLetter<-substr(x,start,start)##get new letter
		if(curLetter!=ignore)
		{
			letterCount<-letterCount+1##conditional increment of letter count	
		}
		res<-paste(res,curLetter,sep="")##push into new peptide sequence
		
		start<-start+1##step forward the pointer after reading one letter
		
		if(letterCount==len)##stop reading letter when letter count reach the pep length
			return(res)
	}
	
}



##read multi_aligned clade sequences from MUSCEL
clade_seqs<-readFASTA(file="~/Programs/git/HIV.db/inst/extdata/MuscleLANL_239_E660.fasta")
cladeNames<-unlist(lapply(clade_seqs,"[[",1))
clade_seqs<-unlist(lapply(clade_seqs,"[[",2))
names(clade_seqs)<-cladeNames
clade_names<-c("mac239", "E660")
mac239_seq<-clade_seqs["mac239_LANL"]

##compute for each clade
for(curCladeName in clade_names)
{
	print(curCladeName)
	
	curSeq<-clade_seqs[curCladeName]#get aligned clade seuqence
	
	#get ordered peptide set from previous step
	curClade<-res[[curCladeName]]
	curPepClade<-curClade$peptide
	nPep<-length(curPepClade)##totoal number of peptides in current clade
	
	#init the gap ranges
#	ir_gap<-vector(mode="list",nPep)
	NoGapRegions<-vector(mode="list",nPep)
	##init standard HXB2 position
	mac239_start<-vector(mode="integer",nPep)	
	mac239_end<-vector(mode="integer",nPep)
	##init aligned global position
	align_start<-vector(mode="integer",nPep)	
	align_end<-vector(mode="integer",nPep)
	newPeps<-vector(mode="character",nPep)
	
	for(i in 1:nPep)#fetch pepetide sequentially from peptide set
	{
		curPep<-curPepClade[i]#get pme peptide
		#compute the new start position
		if(i>1)
		{
			nStep<-curClade$start[i]-curClade$start[i-1]
			preSeq<-substrEx(curSeq,align_start[i-1],nStep)
			nStep_new<-nchar(preSeq)#adjust steps by possible gaps
			align_start[i]<-align_start[i-1]+nStep_new##get new start position for current peptide
		}else
		{
			align_start[i]<-1
		}
		
		print(align_start[i])
		
		##compute the new end position by including gaps
		newPeps[i]<-substrEx(curSeq,align_start[i],nchar(curPep))
		align_end[i]<-align_start[i]+nchar(newPeps[i])-1 
		
		###################################################
		#map to HXB2 standard position
		###################################################
		
		##get HXB2 start by removing gaps
		preSeq<-gsub("-","",substr(mac239_seq,1,align_start[i]-1))
		mac239_start[i]<-nchar(preSeq)+1 
		#get HXB2 end by adding the length of new pep
		mac239_end[i]<-mac239_start[i]+nchar(newPeps[i])-1 
		
		##adjusted the end position by remove the insertions in HXB2 sequence (using gapped ranges)
		HXB2_pep<-substr(mac239_seq,align_start[i],align_end[i])#get peptide sequence from HXB2 aligned sequence
		gapRes<-matchPattern(pattern="-",HXB2_pep)##search for gaps
		mac239_end[i]<-mac239_end[i]-length(gapRes)
		
		
		##save the relative gap positions for for visualization purpose 
		ir_hide<-reduce(ranges(gapRes))
		NoGapRegions[[i]]<-setdiff(IRanges(1,nchar(newPeps[i])),ir_hide)
		names(NoGapRegions)[[i]]<-newPeps[i]
		
#		ir_gap[[i]]<-ir_hide
#		names(ir_gap)[[i]]<-newPeps[i]
		
	}
	res[[curCladeName]]$HXB2_positions<-IRanges(start=mac239_start,end=mac239_end,names=newPeps)	
	res[[curCladeName]]$aligned_positions<-IRanges(start=align_start,end=align_end,names=newPeps)
#	res[[curCladeName]]$gaps<-ir_gap
	res[[curCladeName]]$HXB2GapFilters<-as(IRangesList(NoGapRegions), "GappedRanges")
	
}

######################################################################
##4.merge the new positions from 2 clades
###################################################################
mac239_ir<-IRanges()
ir_mac239GapFilters<-as(IRangesList(),"GappedRanges")
#gap_ir<-NULL
clade<-c()
peptideNb<-c()
aligned<-c()
for(i in 1:nrow(JPTO))
{
	curRow<-JPTO[i,]
	curPep<-as.character(curRow$ID)
#	curCladeName<-clade_names[!is.na(curRow[,clade_names])][1] #list of the clades for this sequence
	curCladeName<-tail(strsplit(as.character(curRow$Name), split="-")[[1]], 1)
	peptideNb<-c(peptideNb,strsplit(as.character(curRow$Name), split="_")[[1]][1])

	if(curCladeName=="E660239")
		{
			curCladeName<-c("mac239", "E660")
			ind<-which(res[["mac239"]]$peptide==curPep)
			mac239_ir<-c(mac239_ir,res[["mac239"]]$HXB2_positions[ind])
			aligned<-c(aligned, res[["mac239"]]$peptide[ind])
			ir_mac239GapFilters<-c(ir_mac239GapFilters,res[["mac239"]]$HXB2GapFilters[ind])
			clade<-c(clade, "mac239,E660")
		}
	else if(curCladeName=="239")
	{
		curCladeName<-"mac239"
		ind<-which(res[[curCladeName]]$peptide==curPep)
		aligned<-c(aligned, res[["mac239"]]$peptide[ind])
		mac239_ir<-c(mac239_ir,res[[curCladeName]]$HXB2_positions[ind])
		ir_mac239GapFilters<-c(ir_mac239GapFilters,res[[curCladeName]]$HXB2GapFilters[ind])
		clade<-c(clade, curCladeName)
	}
	else
	{
		clade<-c(clade, curCladeName)
		ind<-which(res[[curCladeName]]$peptide==curPep)
		curSeq<-res$mac239$consensus
		aligned<-c(aligned, substr(curSeq, start(res$E660$HXB2_positions[ind]), end(res$E660$HXB2_positions[ind])))
		print(c(start(res$E660$HXB2_positions[ind]), end(res$E660$HXB2_positions[ind])))
		mac239_ir<-c(mac239_ir,res[[curCladeName]]$HXB2_positions[ind])
		ir_mac239GapFilters<-c(ir_mac239GapFilters,res[[curCladeName]]$HXB2GapFilters[ind])
	}
#	ind<-which(res[[curCladeName]]$peptide==curPep)
#	mac239_ir<-c(mac239_ir,res[[curCladeName]]$HXB2_positions[ind])
##	gap_ir<-c(gap_ir,res[[curCladeName]]$gaps[ind])
#	ir_mac239GapFilters<-c(ir_mac239GapFilters,res[[curCladeName]]$HXB2GapFilters[ind])
}


######################################################################
##5.construct rangeddata to with the aligned sequence and gapped ranges
###################################################################
#replace the names of range with original pep sequence
names(mac239_ir)<-as.character(JPTO$ID)
#clade<-!is.na(JPTO[,clade_names]) #Instead of the previous clade matrix, we want a simple column indicating the clades

pep_mac239<-RangedData(mac239_ir,aligned,#aligned_Seq,ir_mac239GapFilters,
		peptideNb,clade,space="gp160")
colnames(pep_mac239)[2]<-"gapFilter"
colnames(pep_mac239)[1]<-"aligned"


### get the z-score for each peptide
mat<-.zSum(rownames(pep_mac239[1,]))
for(i in 2:nrow(pep_mac239))
{
	mat<-rbind(mat,.zSum(rownames(pep_mac239[i,])))
}
pep_mac239<-RangedData(ranges(pep_mac239),values(pep_mac239),mat)


######################################################################
##6.compute the trimmed sequences based on aligned sequence and gapped ranges
###################################################################

##trim the sequence by removing the gap
trim_seq<-vector("character",nrow(pep_mac239))
for(i in 1:nrow(pep_mac239))
{
#	browser()
	curSeq<-pep_mac239[["aligned"]][i]
	
	##get the gap boundaries
	curGapFilter<-pep_mac239[["gapFilter"]][i][[1]]
	gaps<-setdiff(IRanges(1,nchar(curSeq)),curGapFilter)
	indToChange<-c(start(gaps)-1,end(gaps)+1)
	
	##flag the start and end of gaps by lower case
	curSeq<-strsplit(curSeq,"")[[1]]
	curSeq[indToChange]<-tolower(curSeq[indToChange])
	
	##trim the flagged sequence
	trim_seq[i]<-paste(seqselect(curSeq,curGapFilter),collapse="")			
}

pep_mac239<-RangedData(ranges(pep_mac239),aligned=pep_mac239[["aligned"]],trimmed=trim_seq,clade)

##add two more columns
anno<-JPTO[,c("Annotation","SEQNO")]
pep_mac239<-RangedData(ranges(pep_mac239),values(pep_mac239),anno)
