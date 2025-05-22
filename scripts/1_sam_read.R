#' Reads SAM file creates a \code{data.frame}
#'
#' This funtion reads SAM file reverts header to columns and outputs
#' a \code{"data.frame"}.
#'
#' @title Read SAM files.
#' @param myfile the name of the fasta file which the data are to be read from.
#' @return a \code{data.frame}.
#' @author Grischa Toedt, Vladislava Milchevskaya \email{milchv@gmail.com}
#' @examples 
#' sam_file <- 
#'  dir(system.file("extdata",package="pdProbeRemap"), 
#'                   pattern="example_drosophila.Aligned.out.sam",
#'                   full.names=TRUE)
#' SAM_list <- read.sam(sam_file)
#' SAM_list$x
#' SAM_list$header
#' @export
###############################################################################
#
# read a sam file and extract the header vars
#
# Author: toedt (but i've changed cat for message, otherwise it would not shut up! Vlada)
###############################################################################
read.sam<- function(myfile){
  message("reading file...")
  x = readLines(myfile)
  message("fetching header...")
  headerpos = grep("^@",x)
  header = x[headerpos]
  message("converting header...")
  header = list("HD" = lapply(gsub("^@HD\t","",header[grep("^@HD",header)]),function(x)strsplit(x,"\t")[[1]]),
                "SQ" = lapply(gsub("^@SQ\t","",header[grep("^@SQ",header)]),function(x)strsplit(x,"\t")[[1]]),
                "RG" = lapply(gsub("^@RG\t","",header[grep("^@RG",header)]),function(x)strsplit(x,"\t")[[1]]),
                "PG" = lapply(gsub("^@PG\t","",header[grep("^@PG",header)]),function(x)strsplit(x,"\t")[[1]]),
                "CO" = gsub("^@CO\t","",header[grep("^@CO",header)])
  )
  message("fetching data...")
  x = lapply(x[-headerpos],line2vals)
  message("convert to data.frame...")
  x = matrix(unlist(x),nrow=length(x),ncol=length(x[[1]]),byrow=TRUE,dimnames=list(NULL,c("QNAME","FLAG","RNAME","POS","MAPQ","CIGAR","RNEXT","PNEXT","TLEN","SEQ","QUAL","ATTR")))
  x = as.data.frame(x,stringsAsFactors=FALSE)
  x$FLAG<- as.integer(x$FLAG)
  x$POS<- as.integer(x$POS)
  x$MAPQ<- as.integer(x$MAPQ)
  x$PNEXT<- as.integer(x$PNEXT)
  x$TLEN<- as.integer(x$TLEN)
  message("done!\n")
  return(list("header"=header,"x" = x))
}
# debug session:
# you could do: take the first 11 fields, the rest are flags
#line2vals<- function(x){
#	c(strsplit(x,"\t",fixed=TRUE)[[1]][1:11],gsub("^([^\t]*[\t]){11}","",x))
#}
# but this is faster faster!
line2vals<- function(x,nfields=11){
  x = strsplit(x,"\t",fixed=TRUE)[[1]]
  c(x[1:nfields],paste(x[-c(1:nfields)],collapse="\t"))
}
# debug session:
#header2mat<- function(x,tag=c("SQ","RG","PG")){
#	tag=match.arg(tag)
#	nfields = c("SQ"=5,"RG"=3,"PG"=3)
#	lapply(gsub(paste("^@",tag,"\t",sep=""),"",x),line2vals,nfields[tag])
#}