library(gdata)
library(Biostrings)
 
JPTdesign<-read.xls("~/Programs/git/HIV.db/inst/extdata/JPT_Annotation_gp160_slide_18Mar11.xls")
#########################################################
##1.pair wise alignment of consecutive peptide sequence to construct entire clade sequence
#######################################################
res<-vector("list",7)
names(res)<-c("M","A","B","C","D","CRF01","CRF02")
for(c in 1:7)
{
  clades<-as.character(JPTdesign[,6+c])
  pepClade<-as.character(JPTdesign$Name[!is.na(clades)])
  start<-rep(0,length(pepClade))
  end<-rep(0,length(pepClade))
  tmpSeq<-do.call(function(x,...){paste(x,sep="",...)},as.list(rep("-",860)))
  start[1]<-1
  end[1]<-start[1]+nchar(pepClade[1])-1
  substr(tmpSeq,start[1],end[1])<-pepClade[1]
  for(i in 2:length(pepClade))
  {
    tmp<-pairwiseAlignment(pepClade[i],pepClade[i-1],type="overlap")
    start[i]<-start[i-1]+start(tmp@subject@range)-1
    end[i]<-start[i]+nchar(pepClade[i])-1
    substr(tmpSeq,start[i],end[i])<-pepClade[i]
  }
  ### I remove the unecessary characters at the end
  res[[c]]<-list(peptide=pepClade,start=start,end=end,consensus=substr(tmpSeq,start[1],end[length(pepClade)]))
}

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
clade_seqs<-readFASTA(file="~/Programs/git/HIV.db/inst/extdata/MuscleMultipleAlignment.fasta")
cladeNames<-unlist(lapply(clade_seqs,"[[",1))
clade_seqs<-unlist(lapply(clade_seqs,"[[",2))
names(clade_seqs)<-cladeNames
clade_names<-c("M","A","B","C","D","CRF01","CRF02")
HXB2_seq<-clade_seqs["HXB2"]

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
	HXB2_start<-vector(mode="integer",nPep)	
	HXB2_end<-vector(mode="integer",nPep)
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
		preSeq<-gsub("-","",substr(HXB2_seq,1,align_start[i]-1))
		HXB2_start[i]<-nchar(preSeq)+1 
		#get HXB2 end by adding the length of new pep
		HXB2_end[i]<-HXB2_start[i]+nchar(newPeps[i])-1 

		##adjusted the end position by remove the insertions in HXB2 sequence (using gapped ranges)
		HXB2_pep<-substr(HXB2_seq,align_start[i],align_end[i])#get peptide sequence from HXB2 aligned sequence
		gapRes<-matchPattern(pattern="-",HXB2_pep)##search for gaps
		HXB2_end[i]<-HXB2_end[i]-length(gapRes)
		 
		
		##save the relative gap positions for for visualization purpose 
		ir_hide<-reduce(ranges(gapRes))
		NoGapRegions[[i]]<-setdiff(IRanges(1,nchar(newPeps[i])),ir_hide)
		names(NoGapRegions)[[i]]<-newPeps[i]
		
#		ir_gap[[i]]<-ir_hide
#		names(ir_gap)[[i]]<-newPeps[i]
		
	}
	res[[curCladeName]]$HXB2_positions<-IRanges(start=HXB2_start,end=HXB2_end,names=newPeps)	
	res[[curCladeName]]$aligned_positions<-IRanges(start=align_start,end=align_end,names=newPeps)
#	res[[curCladeName]]$gaps<-ir_gap
	res[[curCladeName]]$HXB2GapFilters<-as(IRangesList(NoGapRegions), "GappedRanges")
			
}

######################################################################
##4.merge the new positions from 7 clades
###################################################################

##merge the new positions from 7 clades

hxb2_ir<-IRanges()
ir_HXB2GapFilters<-as(IRangesList(),"GappedRanges")
#gap_ir<-NULL
for(i in 1:nrow(JPTdesign))
{
	curRow<-JPTdesign[i,]
	curPep<-as.character(curRow$Name)
	curCladeName<-clade_names[!is.na(curRow[,clade_names])][1]
	ind<-which(res[[curCladeName]]$peptide==curPep)
	hxb2_ir<-c(hxb2_ir,res[[curCladeName]]$HXB2_positions[ind])
#	gap_ir<-c(gap_ir,res[[curCladeName]]$gaps[ind])
	ir_HXB2GapFilters<-c(ir_HXB2GapFilters,res[[curCladeName]]$HXB2GapFilters[ind])
}


######################################################################
##5.construct rangeddata to with the aligned sequence and gapped ranges
###################################################################
#replace the names of range with original pep sequence
names(hxb2_ir)<-as.character(JPTdesign$Name)
clade<-!is.na(JPTdesign[,clade_names])

pep_hxb2<-RangedData(hxb2_ir,aligned_Seq,ir_HXB2GapFilters,clade,space="gp160")
colnames(pep_hxb2)[2]<-"gapFilter"
colnames(pep_hxb2)[1]<-"aligned"



######################################################################
##6.compute the trimmed sequences based on aligned sequence and gapped ranges
###################################################################

##trim the sequence by removing the gap
trim_seq<-vector("character",nrow(pep_hxb2))
for(i in 1:nrow(pep_hxb2))
{
#	browser()
	curSeq<-pep_hxb2[["aligned"]][i]
	
	##get the gap boundaries
	curGapFilter<-pep_hxb2[["gapFilter"]][i][[1]]
	gaps<-setdiff(IRanges(1,nchar(curSeq)),curGapFilter)
	indToChange<-c(start(gaps)-1,end(gaps)+1)

	##flag the start and end of gaps by lower case
	curSeq<-strsplit(curSeq,"")[[1]]
	curSeq[indToChange]<-tolower(curSeq[indToChange])
	
	##trim the flagged sequence
	trim_seq[i]<-paste(seqselect(curSeq,curGapFilter),collapse="")			
}

pep_hxb2<-RangedData(ranges(pep_hxb2),aligned=pep_hxb2[["aligned"]],trimmed=trim_seq,clade)

##add two more columns
anno<-JPTdesign[,c("Annotation","SEQNO")]
pep_hxb2<-RangedData(ranges(pep_hxb2),values(pep_hxb2),anno)

##ADD zpep
Zpep<-read.csv("inst/extdata/peptides_zpep.csv",row.names="peptide")
Zpep<-Zpep[row.names(pep_hxb2),]

pep_hxb2<-RangedData(ranges(pep_hxb2),values(pep_hxb2),Zpep)

save(pep_hxb2,file="data/pep_hxb2.rda")




