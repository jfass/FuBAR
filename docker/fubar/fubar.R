#!/usr/bin/env R

library( Rsamtools )
library( ShortRead )
library( rtracklayer )
library( Rsubread )
library( jsonlite )

# Create genome reference, and read in annotation and protein sequences
GRCh38 <- Rsamtools::FaFile( 'ref/GRCh38.fa' )
GRCh38.info <- seqinfo( GRCh38 )
ann <- rtracklayer::import( 'ref/GRCh38.gtf' )
transcripts.aa <- ( 'ref/gencode.v43.pc_translations.fa' )

# Pull in fusion tables
# Arriba:
results.ar <- data.frame()
numlines <- 0
temp <- suppressWarnings( system2( command='wc', args=c( '-l', 'tmp/trna/arriba/fusions.tsv' ), stdout=T, stderr=T ) )
if ( !( grepl('No such file',temp[1]) ) ) {
  numlines <- as.numeric( sub( " .*", '', temp[1] ) )
}
temp <- data.frame()
if ( numlines > 1 ) {
  temp <- read.csv( file='tmp/trna/arriba/fusions.tsv', sep="\t",
                    header=T, comment.char='', quote='' )
  # fix some column names
  colnames( temp )[c(1,3,4)] <- c('gene1','strand1.gene.fusion','strand2.gene.fusion')
  # create results data frame; fix all "6" --> "chr6" after all tools' results are loaded
  results.ar <- data.frame( tool='Arriba',
                            gene1.name=temp$gene1,
                            gene1.chr=sub(":.*$",'',temp$breakpoint1),
                            gene1.bpt=sub("^.*:",'',temp$breakpoint1),
                            gene1.std=sub("^.*\\/",'',temp$strand1.gene.fusion),
                            gene2.name=temp$gene2,
                            gene2.chr=sub(":.*$",'',temp$breakpoint2),
                            gene2.bpt=sub("^.*:",'',temp$breakpoint2),
                            gene2.std=sub("^.*\\/",'',temp$strand2.gene.fusion),
                            reads=temp$read_identifiers )
  results.ar$gene1.chr[ which( results.ar$gene1.chr=='MT' ) ] <- 'M'
  results.ar$gene1.chr <- paste( 'chr', results.ar$gene1.chr, sep='' )
  results.ar$gene2.chr[ which( results.ar$gene2.chr=='MT' ) ] <- 'M'
  results.ar$gene2.chr <- paste( 'chr', results.ar$gene2.chr, sep='' )
}
# EasyFuse?:
results.ef <- data.frame( tool=character(),
                          gene1.name=character(),
                          gene1.chr=character(),
                          gene1.bpt=character(),
                          gene1.std=character(),
                          gene2.name=character(),
                          gene2.chr=character(),
                          gene2.bpt=character(),
                          gene2.std=character(),
                          reads=character() )  # need to get working after docker refactor
# FusionCatcher:
results.fc <- data.frame()
numlines <- 0
temp <- suppressWarnings( system2( command='wc', args=c( '-l', 'tmp/trna/fusioncatcher/final-list_candidate-fusion-genes.txt' ),
                                   stdout=T, stderr=T ) )
if ( !( grepl('No such file',temp[1]) ) ) {
  numlines <- as.numeric( sub( " .*", '', temp[1] ) )
}
temp <- data.frame()
if ( numlines > 1 ) {
  temp <- read.csv( file='tmp/trna/fusioncatcher/final-list_candidate-fusion-genes.txt', sep="\t",
                    header=T, comment.char='', quote='' )
  # fix some column names
  colnames( temp ) <- sub( "\\.$", '', colnames( temp ) )
  # create results data frame; fix all "6" --> "chr6" after all tools' results are loaded
  results.fc <- data.frame( tool='FusionCatcher',
                            gene1.name=temp$Gene_1_symbol.5end_fusion_partner,
                            gene1.chr=sub(":.*$",'',temp$Fusion_point_for_gene_1.5end_fusion_partner),
                            gene1.bpt=sub(":.*$",'',sub("^.*?:",'',temp$Fusion_point_for_gene_1.5end_fusion_partner)),
                            gene1.std=sub("^.*:",'',temp$Fusion_point_for_gene_1.5end_fusion_partner),
                            gene2.name=temp$Gene_2_symbol.3end_fusion_partner,
                            gene2.chr=sub(":.*$",'',temp$Fusion_point_for_gene_2.3end_fusion_partner),
                            gene2.bpt=sub(":.*$",'',sub("^.*?:",'',temp$Fusion_point_for_gene_2.3end_fusion_partner)),
                            gene2.std=sub("^.*:",'',temp$Fusion_point_for_gene_2.3end_fusion_partner),
                            reads='' )  # reads are renamed, in own files, in tmp/fc/
  results.fc$gene1.chr[ which( results.fc$gene1.chr=='MT' ) ] <- 'M'
  results.fc$gene1.chr <- paste( 'chr', results.fc$gene1.chr, sep='' )
  results.fc$gene2.chr[ which( results.fc$gene2.chr=='MT' ) ] <- 'M'
  results.fc$gene2.chr <- paste( 'chr', results.fc$gene2.chr, sep='' )
}
# MapSplice:
results.ms <- data.frame()
numlines <- 0
temp <- suppressWarnings( system2( command='wc', args=c( '-l', 'tmp/trna/mapsplice/fusions_well_annotated.txt' ),
                                   stdout=T, stderr=T ) )
if ( !( grepl('No such file',temp[1]) ) ) {
  numlines <- as.numeric( sub( " .*", '', temp[1] ) )
}
if ( numlines > 0 ) {
  temp <- read.csv( file='tmp/trna/mapsplice/fusions_well_annotated.txt', sep='\t',
                    header=F, comment.char='', quote='' )
  # add column names (!)
  colnames( temp ) <- c( 'chrom', 'doner_end', 'acceptor_start', 'id', 'coverage', 'strand', 'rgb', 
                         'block_count', 'block_size', 'block_distance', 'entropy', 'flank_case', 'flank_string',
                         'min_mismatch', 'max_mismatch', 'ave_mismatch',
                         'max_min_suffix', 'max_min_prefix', 'min_anchor_difference',
                         'unique_read_count', 'multi_read_count', 'paired_read_count',
                         'left_paired_read_count', 'right_paired_read_count',
                         'multiple_paired_read_count', 'unique_paired_read_count',
                         'single_read_count', 'encompassing_readpair_count',
                         'doner_start', 'acceptor_end', 'doner_isoforms', 'acceptor_isoforms',
                         'dus_obsolete', 'aus_obsolete', 'duKSs_obsolete', 'auKSs_obsolete',
                         'minimal_doner_isoform_length', 'maximal_doner_isoform_length',
                         'minimal_acceptor_isoform_length', 'maximal_acceptor_isoform_length',
                         'paired_reads_entropy', 'mismatch_per_bp', 'anchor_score',
                         'max_doner_fragment', 'max_acceptor_fragment',
                         'max_cur_fragment', 'min_cur_fragment', 'ave_cur_fragment',
                         'doner_encompass_unique', 'doner_encompass_multiple',
                         'acceptor_encompass_unique', 'acceptor_encompass_multiple',
                         'doner_match_to_normal', 'acceptor_match_to_normal',
                         'doner_seq', 'acceptor_seq',
                         'match_gene_strand', 'annotated_type', 'fusion_type', 'gene_strand',
                         'annotated_gene_donor', 'annotated_gene_acceptor' )
  # create results data frame
  results.ms <- data.frame( tool='MapSplice',
                            gene1.name=sub(",.*$",'',temp$annotated_gene_donor),
                            gene1.chr=sub("~.*$",'',temp$chrom),
                            gene1.bpt=temp$doner_end,
                            gene1.std=sub( ".$", '', temp$strand ),
                            gene2.name=sub(",.*$",'',temp$annotated_gene_acceptor),
                            gene2.chr=sub("^.*~",'',temp$chrom),
                            gene2.bpt=temp$acceptor_start,
                            gene2.std=sub( "^.", '', temp$strand ),
                            reads='' )  # need to grep out of alignment.sam/bam
  temp <- suppressWarnings( system2( command='grep', args='-F "ZF:Z:FUS_" tmp/trna/mapsplice/alignments.sam',
                            stdout=T, stderr=T ) )
  for (i in 1:dim(results.ms)[1]) {
    fus <- paste( 'FUS', results.ms$gene1.bpt[i], results.ms$gene2.bpt[i], sep='_' )
    results.ms$reads[ i ] <- paste( unique( sub( "\t.*", '', temp[ grep( fus, temp ) ] ) ), collapse=',' )
  }
}
# STAR-Fusion:
results.sf <- data.frame()
numlines <- 0
temp <- suppressWarnings( system2( command='wc', args=c( '-l', 'tmp/trna/starfusion/star-fusion.fusion_predictions.tsv' ),
                                   stdout=T, stderr=T ) )
if ( !( grepl('No such file',temp[1]) ) ) {
  numlines <- as.numeric( sub( " .*", '', temp[1] ) )
}
temp <- data.frame()
if ( numlines > 1 ) {
  temp <- read.csv( file='tmp/trna/starfusion/star-fusion.fusion_predictions.tsv', sep='\t',
                    header=T, comment.char='', quote='' )
  # fix first column name
  colnames( temp )[1] <- 'FusionName'
  # create results data frame
  results.sf <- data.frame( tool='STAR-Fusion',
                            gene1.name=sub("\\^.*$",'',temp$LeftGene),
                            gene1.chr=sub(":.*$",'',temp$LeftBreakpoint),
                            gene1.bpt=sub(":.*$",'',sub("^.*?:",'',temp$LeftBreakpoint)),
                            gene1.std=sub("^.*:",'',temp$LeftBreakpoint),
                            gene2.name=sub("\\^.*$",'',temp$RightGene),
                            gene2.chr=sub(":.*$",'',temp$RightBreakpoint),
                            gene2.bpt=sub(":.*$",'',sub("^.*?:",'',temp$RightBreakpoint)),
                            gene2.std=sub("^.*:",'',temp$RightBreakpoint),
                            reads=paste(temp$JunctionReads,temp$SpanningFrags,sep=',') )
}
# STAR-SEQR:
results.ss <- data.frame()
numlines <- 0
temp <- suppressWarnings( system2( command='wc',
                                   args=c( '-l', 'tmp/trna/starseqr/sample_STAR-SEQR/sample_STAR-SEQR_candidates.txt' ),
                                   stdout=T, stderr=T ) )
if ( !( grepl('No such file',temp[1]) ) ) {
  numlines <- as.numeric( sub( " .*", '', temp[1] ) )
}
temp <- data.frame()
if ( numlines > 1 ) {
  temp <- read.csv( file='tmp/trna/starseqr/sample_STAR-SEQR/sample_STAR-SEQR_candidates.txt', sep='\t',
                    header=T, comment.char='', quote='' )
  # create results data frame (note STAR-SEQR breakpoints adjusted from zero-based coords)
  results.ss <- data.frame( tool='STAR-SEQR',
                            gene1.name=temp$LEFT_SYMBOL,
                            gene1.chr=sub(":.*$",'',temp$BRKPT_LEFT),
                            gene1.bpt=as.numeric(sub(":.*$",'',sub("^.*?:",'',temp$BRKPT_LEFT)))+1,
                            gene1.std=sub("^.*:",'',temp$BRKPT_LEFT),
                            gene2.name=temp$RIGHT_SYMBOL,
                            gene2.chr=sub(":.*$",'',temp$BRKPT_RIGHT),
                            gene2.bpt=as.numeric(sub(":.*$",'',sub("^.*?:",'',temp$BRKPT_RIGHT)))+1,
                            gene2.std=sub("^.*:",'',temp$BRKPT_RIGHT),
                            reads='' )
  ids.ss <- temp$ID  # internal id consistent with dirnames in tmp/ss/sampl_STAR-SEQR/support/
}
if ( exists( 'ids.ss' ) ) {
  for ( i in 1:length(ids.ss) ) {
    supportregex <- gsub( ':', '_', ids.ss[ i ] )
    supportregex <- gsub( '\\+', 'pos', supportregex )
    supportregex <- gsub( '\\-', 'neg', supportregex )
    dn <- list.files( path='tmp/trna/starseqr/sample_STAR-SEQR/support', pattern=supportregex, include.dirs=T, full.names=T )
    if ( length( dn ) > 0 ) {
      reads.over <- readFastq( paste( dn, '/overhang.fastq', sep='' ) )
      ids.over <- sub( "_.*$", '', as.character( id( reads.over ) ) )  # remove _flag
      reads.span <- readFastq( paste( dn, '/span.fastq', sep='' ) )
      ids.span <- sub( "\\/.*$", '', as.character( id( reads.span ) ) )  # remove /1 or /2
      reads.split <- readFastq( paste( dn, '/split.fastq', sep='' ) )
      ids.split <- sub( "_.*$", '', as.character( id( reads.split ) ) )  # remove _flag
      results.ss$reads[i] <- paste( unique( c( ids.over, ids.span, ids.split ) ), collapse=',' )
    }
  }
}
# "preset" fusions of interest, regardless of evidence or call
results.ps <- data.frame()
numlines <- 0
temp <- suppressWarnings( system2( command='wc',
                                   args=c( '-l', 'in/trna/preset.fusions.tsv' ),
                                   stdout=T, stderr=T ) )
if ( !( grepl('No such file',temp[1]) ) ) {
  numlines <- as.numeric( sub( " .*", '', temp[1] ) )
}
temp <- data.frame()
if ( numlines > 1 ) {
  temp <- read.csv( file='in/trna/preset.fusions.tsv', sep='\t',
                    header=T, comment.char='', quote='' )
  # create results data frame; no read ids because the fusion may not have been called by any tools!
  results.ps <- data.frame( tool='preset',
                            gene1.name=temp$gene1,
                            gene1.chr=temp$chr1,
                            gene1.bpt=temp$bpt1,
                            gene1.std=temp$std1,
                            gene2.name=temp$gene2,
                            gene2.chr=temp$chr2,
                            gene2.bpt=temp$bpt2,
                            gene2.std=temp$std2,
                            reads='' )
  results.ps <- unique( results.ps )
}

# consolidate, track which tools called each fusion
results.all <- rbind( results.ar, results.ef, results.fc,
                      results.ms, results.sf, results.ss,
                      results.ps )
results.all$key <- apply( results.all[,2:9], 1, function(x) {paste(x,collapse='|')} )
uks <- unique( results.all$key )
results <- data.frame( gene1.name=character(),
		       gene1.chr=character(),
		       gene1.bpt=numeric(),
		       gene1.std=character(),
		       gene2.name=character(),
		       gene2.chr=character(),
		       gene2.bpt=numeric(),
		       gene2.std=character(),
		       callers=character(),
		       numcalls=numeric(),
		       reads=character() )
newres <- results
if ( length(uks) > 0 ) {
  for ( k in uks ) {
    #print( k )
    temp <- strsplit( k, split='|', fixed=T )[[1]]
    newres[ 1, ] <- NA
    newres$gene1.name <- temp[ 1 ]
    newres$gene1.chr <- temp[ 2 ]
    newres$gene1.bpt <- as.numeric( temp[ 3 ] )
    newres$gene1.std <- temp[ 4 ]
    newres$gene2.name <- temp[ 5 ]
    newres$gene2.chr <- temp[ 6 ]
    newres$gene2.bpt <- as.numeric( temp[ 7 ] )
    newres$gene2.std <- temp[ 8 ]
    matches <- which( results.all$key == k )
    newres$callers <- paste( results.all$tool[ matches ], collapse=',' )
    newres$numcalls <- length( matches )
    if ( 'preset' %in% results.all$tool[matches] ) { newres$numcalls <- newres$numcalls - 1 }  # include 'preset' in caller list, but don't count it
    temp.all.reads <- paste( results.all$reads[ matches ], collapse=',' )
    if ( length( temp.all.reads ) > 0 ) {
      temp.all.reads <- unlist( strsplit( temp.all.reads, ',', fixed=T ) )
      # only include read ids that conform to expected Illumina colon-separated format
      temp.all.reads <- temp.all.reads[ grep( "^.*:.*:.*:.*:.*:.*:.*$",
                                              temp.all.reads ) ]
      newres$reads <- paste( temp.all.reads, collapse=',' )
    }
    results <- rbind( results, newres )
  }
  results$support_count <- as.numeric( NA )
}

# load subsets of read sets attributed to fusions (by callers) or aligned adjacent to breakpoints
# ... build regions lists for selecting reads for realignment, bracketing breakpoints by 60nt to help inform phasing of nearby variants:
regions1.gr <- GRanges( seqnames = results$gene1.chr,
                        ranges=IRanges( start = results$gene1.bpt-60,
                                        end = results$gene1.bpt+60 ) )
for ( i in 1:length( regions1.gr ) ) {
  if ( start(regions1.gr)[i] < 1 ) { start( regions1.gr)[i] <- 1 }  # too far left; peg to 1
  if ( end( regions1.gr)[i] > seqlengths( GRCh38.info )[ which( seqnames( GRCh38.info ) == seqnames( regions1.gr )[i] ) ] ) {
    end( regions1.gr )[i] <- seqlengths( GRCh38.info )[ which( seqnames( GRCh38.info ) == seqnames( regions1.gr )[i] ) ]
  }  # too far right; peg to chromosome length
}
regions2.gr <- GRanges( seqnames = results$gene2.chr,
                        ranges=IRanges( start = results$gene2.bpt-60,
                                        end = results$gene2.bpt+60 ) )
for ( i in 1:length( regions2.gr ) ) {
  if ( start(regions2.gr)[i] < 1 ) { start( regions2.gr)[i] <- 1 }  # too far left; peg to 1
  if ( end( regions2.gr)[i] > seqlengths( GRCh38.info )[ which( seqnames( GRCh38.info ) == seqnames( regions2.gr )[i] ) ] ) {
    end( regions2.gr )[i] <- seqlengths( GRCh38.info )[ which( seqnames( GRCh38.info ) == seqnames( regions2.gr )[i] ) ]
  }  # too far right; peg to chromosome length
}
# ... collapse to sorted list:
seqlevels( regions1.gr ) <- paste( 'chr', c( as.character( 1:22), 'X', 'Y' ), sep='' )
seqlevels( regions2.gr ) <- paste( 'chr', c( as.character( 1:22), 'X', 'Y' ), sep='' )
reduced.gr <- sort( reduce( c( regions1.gr, regions2.gr ) ) )  # sort, uniquify, and collapse overlapping ranges
reduced.list <- paste( seqnames(reduced.gr), ':', start(reduced.gr), '-', end(reduced.gr), sep='', collapse=' ' )  # list for 'samtools view' command
# put together read id sets for rnaseq and whole exome samples:
# ... trna:
temp <- suppressWarnings( system2( command='samtools',
                                   args=c( 'view',
                                           '-@ 6',
                                           'tmp/trna/recal.bam',
                                           reduced.list,
                                           '| cut -f1 | sort -u > tmp/trna/fusionqc/name.lst' ),
                                   stdout=T, stderr=T ) )  # create list of read ids from main trna bam
temp <- suppressWarnings( system2( command='samtools',
                                   args=c( 'view',
                                           '-@ 6',
                                           'tmp/trna/arriba/Aligned.sortedByCoord.out.bam',
                                           gsub( 'chr', '', reduced.list ),
                                           '| cut -f1 | sort -u >> tmp/trna/fusionqc/name.lst' ),
                                   stdout=T, stderr=T ) )  # add list of read ids from Arriba's bam
temp <- suppressWarnings( system2( command='samtools',
                                   args=c( 'view',
                                           '-@ 6',
                                           'tmp/trna/mapsplice/alignments.bam',
                                           reduced.list,
                                           '| cut -f1 | sort -u >> tmp/trna/fusionqc/name.lst' ),
                                   stdout=T, stderr=T ) )  # add list of read ids from Mapsplice's bam
temp <- suppressWarnings( system2( command='samtools',
                                   args=c( 'view',
                                           '-@ 6',
                                           'tmp/trna/starfusion/Aligned.out.bam',
                                           reduced.list,
                                           '| cut -f1 | sort -u >> tmp/trna/fusionqc/name.lst' ),
                                   stdout=T, stderr=T ) )  # add list of read ids from StarFusion's bam
temp <- suppressWarnings( system2( command='samtools',
                                   args=c( 'view',
                                           '-@ 6',
                                           'tmp/trna/starseqr/sample_STAR-SEQR/sample.Chimeric.out.bam',
                                           reduced.list,
                                           '| cut -f1 | sort -u >> tmp/trna/fusionqc/name.lst' ),
                                   stdout=T, stderr=T ) )  # add list of read ids from Star-SEQR's bam
temp <- suppressWarnings( system2( command='cat',
                                   args=c( 'tmp/trna/fusionqc/name.lst | sort -u > tmp/trna/fusionqc/temp.lst' ),
                                   stdout=T, stderr=T ) )
temp <- suppressWarnings( system2( command='mv',
                                   args=c( 'tmp/trna/fusionqc/temp.lst',
                                           'tmp/trna/fusionqc/name.lst' ),
                                   stdout=T, stderr=T ) )
temp <- suppressWarnings( system2( command='seqtk',
                                   args=c( 'subseq',
                                           'in/trna/r1.fq.gz',
                                           'tmp/trna/fusionqc/name.lst',
                                           '| gzip > tmp/trna/fusionqc/trna.r1.fq.gz' ),
                                   stdout=T, stderr=T ) )  # retrieve read 1
temp <- suppressWarnings( system2( command='seqtk',
                                   args=c( 'subseq',
                                           'in/trna/r2.fq.gz',
                                           'tmp/trna/fusionqc/name.lst',
                                           '| gzip > tmp/trna/fusionqc/trna.r2.fq.gz' ),
                                   stdout=T, stderr=T ) )  # retrieve read 2
# ... twes:
temp <- suppressWarnings( system2( command='samtools',
                                   args=c( 'view',
                                           '-@ 6',
                                           'tmp/twes/recal.bam',
                                           reduced.list,
                                           '| cut -f1 | sort -u > tmp/trna/fusionqc/name.lst' ),
                                   stdout=T, stderr=T ) )  # create list of read ids
temp <- suppressWarnings( system2( command='seqtk',
                                   args=c( 'subseq',
                                           'in/twes/r1.fq.gz',
                                           'tmp/trna/fusionqc/name.lst',
                                           '| gzip > tmp/trna/fusionqc/twes.r1.fq.gz' ),
                                   stdout=T, stderr=T ) )  # retrieve read 1
temp <- suppressWarnings( system2( command='seqtk',
                                   args=c( 'subseq',
                                           'in/twes/r2.fq.gz',
                                           'tmp/trna/fusionqc/name.lst',
                                           '| gzip > tmp/trna/fusionqc/twes.r2.fq.gz' ),
                                   stdout=T, stderr=T ) )  # retrieve read 2
# ... nwes:
temp <- suppressWarnings( system2( command='samtools',
                                   args=c( 'view',
                                           '-@ 6',
                                           'tmp/nwes/recal.bam',
                                           reduced.list,
                                           '| cut -f1 | sort -u > tmp/trna/fusionqc/name.lst' ),
                                   stdout=T, stderr=T ) )  # create list of read ids
temp <- suppressWarnings( system2( command='seqtk',
                                   args=c( 'subseq',
                                           'in/nwes/r1.fq.gz',
                                           'tmp/trna/fusionqc/name.lst',
                                           '| gzip > tmp/trna/fusionqc/nwes.r1.fq.gz' ),
                                   stdout=T, stderr=T ) )  # retrieve read 1
temp <- suppressWarnings( system2( command='seqtk',
                                   args=c( 'subseq',
                                           'in/nwes/r2.fq.gz',
                                           'tmp/trna/fusionqc/name.lst',
                                           '| gzip > tmp/trna/fusionqc/nwes.r2.fq.gz' ),
                                   stdout=T, stderr=T ) )  # retrieve read 2
# pull in (reduced set of) trna reads, to extract and remap to fusion qc references
trna.R1 <- readFastq( 'tmp/trna/fusionqc/trna.r1.fq.gz' )
trna.R2 <- readFastq( 'tmp/trna/fusionqc/trna.r2.fq.gz' )
trna.ids <- sub( "\\ .*$", '', as.character( id( trna.R1 ) ) )
# pull in (reduced set of) nwes reads, to extract and remap to fusion qc references
nwes.R1 <- readFastq( 'tmp/trna/fusionqc/nwes.r1.fq.gz' )
nwes.R2 <- readFastq( 'tmp/trna/fusionqc/nwes.r2.fq.gz' )
nwes.ids <- sub( "\\ .*$", '', as.character( id( nwes.R1 ) ) )
# pull in (reduced set of) twes reads, to extract and remap to fusion qc references
twes.R1 <- readFastq( 'tmp/trna/fusionqc/twes.r1.fq.gz' )
twes.R2 <- readFastq( 'tmp/trna/fusionqc/twes.r2.fq.gz' )
twes.ids <- sub( "\\ .*$", '', as.character( id( twes.R1 ) ) )

# loop through fusions, creating qc ref and remapping reads
for (i in 1:dim(results)[1]) {
  # retrieve padded gene1 genomic sequence + annotation, and adjust
  gene1.name <- results$gene1.name[i]
  gene1.chr <- results$gene1.chr[i]
  gene1.std <- results$gene1.std[i]  # fusion strand, not gene strand
  gene1.bpt <- as.numeric( results$gene1.bpt[i] )
  gene1.ann <- ann[ which( ann$gene_name == gene1.name &
			   as.character(seqnames(ann)) == gene1.chr &
			   as.character(strand(ann)) == gene1.std ) ]
  if ( length(gene1.ann) > 0 ) {
    gene1.code <- names( sort( table( strand( gene1.ann ) ), decreasing=T )[ 1 ] )  # coding strand
    gene1.left <- min( start( gene1.ann ), gene1.bpt ) - 5000  # note bpt could be outside gene boundaries
    if ( gene1.left < 1 ) { gene1.left <- 1 }  # too far left; peg to 1
    gene1.right <- max( end( gene1.ann ), gene1.bpt ) + 5000
    if ( gene1.right > seqlengths( GRCh38.info )[ which( seqnames( GRCh38.info ) == gene1.chr ) ] ) {
      gene1.right <- seqlengths( GRCh38.info )[ which( seqnames( GRCh38.info ) == gene1.chr ) ]
    }  # too far right; peg to ref sequence's length
    gene1.gr <- GRanges( seqnames=gene1.chr, strand='+',
                         ranges=IRanges( start=gene1.left, end=gene1.right ) )
    gene1.seq <- Rsamtools::scanFa( GRCh38, gene1.gr )
    names( gene1.seq ) <- paste( 'gene1', gene1.name, sep='_' )
    offset <- gene1.left - 1
    temp.uniqname <- paste( 'gene1', gene1.name, sep='_' )  # for when gene1=gene2
    temp.gr <- GRanges( seqnames=rep(temp.uniqname,length(gene1.ann)),
                        strand=strand(gene1.ann),
                        ranges=IRanges(start=start(gene1.ann)-offset,end=end(gene1.ann)-offset) )
    mcols( temp.gr ) <- mcols( gene1.ann )
    # identify fusion partner fragment
    fuse1.gr <- GRanges( seqnames=temp.uniqname, strand='*',
                         ranges=IRanges(start=1,end=gene1.right-offset) )
    fuse1.gr$gene_name <- 'fusion_partner_1'
    if ( gene1.std=='+' ) {
      strand( fuse1.gr ) <- '+'
      end( fuse1.gr ) <- gene1.bpt - offset
    } else if ( gene1.std=='-' ) {
      strand( fuse1.gr ) <- '-'
      start( fuse1.gr ) <- gene1.bpt - offset
    }
    gene1.ann.orig <- gene1.ann  # hand off annotation in genomic coordinates
    gene1.ann <- temp.gr  # new annotation in local coordinates
  }
  # retrieve padded gene2 genomic sequence + annotation, and adjust
  gene2.name <- results$gene2.name[i]
  gene2.chr <- results$gene2.chr[i]
  gene2.std <- results$gene2.std[i]  # fusion strand, not gene strand
  gene2.bpt <- as.numeric( results$gene2.bpt[i] )
  gene2.ann <- ann[ which( ann$gene_name == gene2.name &
			   as.character(seqnames(ann)) == gene2.chr &
			   as.character(strand(ann)) == gene2.std ) ]
  if ( length(gene2.ann) > 0 ) {
    gene2.code <- names( sort( table( strand( gene2.ann ) ), decreasing=T )[ 1 ] )  # coding strand
    gene2.left <- min( start( gene2.ann ), gene2.bpt ) - 5000  # note bpt could be outside gene boundaries
    if ( gene2.left < 1 ) { gene2.left <- 1 }  # too far left; peg to 1
    gene2.right <- max( end( gene2.ann ), gene2.bpt ) + 5000
    if ( gene2.right > seqlengths( GRCh38.info )[ which( seqnames( GRCh38.info ) == gene2.chr ) ] ) {
      gene2.right <- seqlengths( GRCh38.info )[ which( seqnames( GRCh38.info ) == gene2.chr ) ]
    }  # too far right; peg to ref sequence's length
    gene2.gr <- GRanges( seqnames=gene2.chr, strand='+',
                         ranges=IRanges( start=gene2.left, end=gene2.right ) )
    gene2.seq <- Rsamtools::scanFa( GRCh38, gene2.gr )
    names( gene2.seq ) <- paste( 'gene2', gene2.name, sep='_' )
    offset <- gene2.left - 1
    temp.uniqname <- paste( 'gene2', gene2.name, sep='_' )  # for when gene1=gene2
    temp.gr <- GRanges( seqnames=rep(temp.uniqname,length(gene2.ann)),
                        strand=strand(gene2.ann),
                        ranges=IRanges(start=start(gene2.ann)-offset,end=end(gene2.ann)-offset) )
    mcols( temp.gr ) <- mcols( gene2.ann )
    fuse2.gr <- GRanges( seqnames=temp.uniqname, strand='*',
                         ranges=IRanges(start=1,end=gene2.right-offset) )
    fuse2.gr$gene_name <- 'fusion_partner_2'
    if ( gene2.std=='+' ) {
      strand( fuse2.gr ) <- '+'
      start( fuse2.gr ) <- gene2.bpt - offset
    } else if ( gene2.std=='-' ) {
      strand( fuse2.gr ) <- '-'
      end( fuse2.gr ) <- gene2.bpt - offset
    }
    gene2.ann.orig <- gene2.ann  # hand off annotation in genomic coordinates
    gene2.ann <- temp.gr  # new annotation in local coordinates
  }
  # construct fusion and its chimeric annotation, depending on strands
  if ( ( length( gene1.ann ) > 0 ) & ( length( gene2.ann ) > 0 ) ) {
    if (gene1.std == '+') {
      fusion1.gr <- GRanges( seqnames=gene1.chr, strand='+',
                             ranges=IRanges( start=gene1.left, end=gene1.bpt ) )
      fusion1.seq <- Rsamtools::scanFa( GRCh38, fusion1.gr )
      fusion1.ann <- gene1.ann.orig[ subjectHits( findOverlaps(fusion1.gr,gene1.ann.orig,ignore.strand=T) ) ]
      drop <- which( start( fusion1.ann ) >= gene1.bpt )  # exclude features starting to right of breakpoint
      if ( length( drop ) > 0 ) { fusion1.ann <- fusion1.ann[ -drop ] }
      end( fusion1.ann )[ which( end( fusion1.ann ) > gene1.bpt ) ] <- gene1.bpt
      offset <- gene1.left - 1
      start( fusion1.ann ) <- start( fusion1.ann ) - offset
      end( fusion1.ann ) <- end( fusion1.ann ) - offset
    } else if (gene1.std == '-') {
      fusion1.gr <- GRanges( seqnames=gene1.chr, strand='+',
                             ranges=IRanges( start=gene1.bpt, end=gene1.right ) )
      fusion1.seq <- Rsamtools::scanFa( GRCh38, fusion1.gr )
      fusion1.seq <- reverseComplement( fusion1.seq )  # GRanges strand doesn't affect scanFa
      fusion1.ann <- gene1.ann.orig[ subjectHits( findOverlaps(fusion1.gr,gene1.ann.orig,ignore.strand=T) ) ]
      drop <-  which( end( fusion1.ann ) <= gene1.bpt )  # exclude features ending to left of breakpoint
      if ( length( drop ) > 0 ) { fusion1.ann <- fusion1.ann[ -drop ] }
      start( fusion1.ann )[ which( start( fusion1.ann ) < gene1.bpt ) ] <- gene1.bpt
      ranges( fusion1.ann ) <- reflect( ranges( fusion1.ann ), bounds=IRanges(start=1,end=gene1.right) )
      if ( all( strand( fusion1.ann ) == '+' ) ) { strand( fusion1.ann ) <- '-' } else
      if ( all( strand( fusion1.ann ) == '-' ) ) { strand( fusion1.ann ) <- '+' }
    }
    if (gene2.std == '+') {
      fusion2.gr <- GRanges( seqnames=gene2.chr, strand='+',
                             ranges=IRanges( start=gene2.bpt, end=gene2.right ) )
      fusion2.seq <- Rsamtools::scanFa( GRCh38, fusion2.gr )
      fusion2.ann <- gene2.ann.orig[ subjectHits( findOverlaps(fusion2.gr,gene2.ann.orig,ignore.strand=T) ) ]
      drop <- which( end( fusion2.ann ) <= gene2.bpt )  # exclude features ending to left of breakpoint
      if ( length( drop ) > 0 ) { fusion2.ann <- fusion2.ann[ -drop ] }
      start( fusion2.ann )[ which( start( fusion2.ann ) < gene2.bpt ) ] <- gene2.bpt
      offset <- gene2.bpt - 1
      start( fusion2.ann ) <- start( fusion2.ann ) - offset + width( fusion1.seq )
      end( fusion2.ann ) <- end( fusion2.ann ) - offset + width( fusion1.seq )
    } else if (gene2.std == '-') {
      fusion2.gr <- GRanges( seqnames=gene2.chr, strand='-',
                             ranges=IRanges( start=gene2.left, end=gene2.bpt ) )
      fusion2.seq <- Rsamtools::scanFa( GRCh38, fusion2.gr )
      fusion2.seq <- reverseComplement( fusion2.seq )  # GRanges strand doesn't affect scanFa
      fusion2.ann <- gene2.ann.orig[ subjectHits( findOverlaps(fusion2.gr,gene2.ann.orig,ignore.strand=T) ) ]
      drop <- which( start( fusion2.ann ) >= gene2.bpt )  # exclude features starting to right of breakpoint
      if ( length( drop ) > 0 ) { fusion2.ann <- fusion2.ann[ -drop ] }
      end( ranges( fusion2.ann ) )[ which( end( ranges( fusion2.ann ) ) > gene2.bpt ) ] <- gene2.bpt
      ranges( fusion2.ann ) <- shift( reflect( ranges( fusion2.ann ), bounds=IRanges(start=1,end=gene2.bpt) ), width( fusion1.seq ) )
      if ( all( strand( fusion2.ann ) == '+' ) ) { strand( fusion2.ann ) <- '-' } else
      if ( all( strand( fusion2.ann ) == '-' ) ) { strand( fusion2.ann ) <- '+' }
    }
    temp.seq <- as.character( paste( fusion1.seq, fusion2.seq, sep='' ) )
    pair.underscores <- paste( gene1.name, gene2.name, sep='__' )
    names( temp.seq ) <- pair.underscores
    fusion.seq <- DNAStringSet( temp.seq )
    temp.ann <- c( fusion1.ann, fusion2.ann )
    # chimeric gene annotation:
    fusion.ann <- GRanges( seqnames=Rle( values=pair.underscores, lengths=length(temp.ann) ),
                           strand=strand( temp.ann ),
                           ranges=ranges( temp.ann ) )
    mcols( fusion.ann ) <- mcols( temp.ann )
    # chimeric fusion partner annotation:
    fused.gr <- GRanges( seqnames=Rle( values=pair.underscores, lengths=2 ),
                         strand=c('+','+'),
                         ranges=IRanges( start=c( 1, seqlengths(fusion1.seq) + 1 ),
                                         end=c( seqlengths(fusion1.seq), seqlengths(fusion1.seq)+seqlengths(fusion2.seq) ) ) )
    fused.gr$gene_name <- c( 'fusion_partner_1', 'fusion_partner_2' )
    # calculate number of supporting trna read pairs, to go into folder naming
    ids.support <- unlist( strsplit( results$reads[ i ], ',' ) )
    trna.ind <- which( trna.ids %in% ids.support )
    results$support_count[ i ] <- length( trna.ind )
    # tack on any reads near breakpoints, bracketing by 60nt to inform neighboring variants that might need phasing
    bpt1.L <- start( regions1.gr )[ i ]
    bpt1.R <- end( regions1.gr )[ i ]
    bpt2.L <- start( regions2.gr )[ i ]
    bpt2.R <- end( regions2.gr )[ i ]
    ids.maybe.ar <- suppressWarnings( system2( command='samtools',
                                               args=c( 'view',
                                                       '-@ 12',
                                                       'tmp/trna/arriba/Aligned.sortedByCoord.out.bam',
                                                       paste(paste(sub('chr','',gene1.chr),bpt1.L,sep=':'),bpt1.R,sep='-'),
                                                       paste(paste(sub('chr','',gene2.chr),bpt2.L,sep=':'),bpt2.R,sep='-'),
                                                       '| cut -f1 | sort -u' ),
                                               stdout=T, stderr=T ) )
    ids.maybe.ms <- suppressWarnings( system2( command='samtools',
                                               args=c( 'view',
                                                       '-@ 12',
                                                       'tmp/trna/mapsplice/alignments.bam',
                                                       paste(paste(gene1.chr,bpt1.L,sep=':'),bpt1.R,sep='-'),
                                                       paste(paste(gene2.chr,bpt2.L,sep=':'),bpt2.R,sep='-'),
                                                       '| cut -f1 | sort -u' ),
                                               stdout=T, stderr=T ) )
    ids.maybe.sf <- suppressWarnings( system2( command='samtools',
                                               args=c( 'view',
                                                       '-@ 12',
                                                       'tmp/trna/starfusion/Aligned.out.bam',
                                                       paste(paste(gene1.chr,bpt1.L,sep=':'),bpt1.R,sep='-'),
                                                       paste(paste(gene2.chr,bpt2.L,sep=':'),bpt2.R,sep='-'),
                                                       '| cut -f1 | sort -u' ),
                                               stdout=T, stderr=T ) )
    ids.maybe.ss <- suppressWarnings( system2( command='samtools',
                                               args=c( 'view',
                                                       '-@ 12',
                                                       'tmp/trna/starseqr/sample_STAR-SEQR/sample.Chimeric.out.bam',
                                                       paste(paste(gene1.chr,bpt1.L,sep=':'),bpt1.R,sep='-'),
                                                       paste(paste(gene2.chr,bpt2.L,sep=':'),bpt2.R,sep='-'),
                                                       '| cut -f1 | sort -u' ),
                                               stdout=T, stderr=T ) )
    ids.maybe <- unique( c( ids.maybe.ar, ids.maybe.ms, ids.maybe.sf, ids.maybe.ss ) )
    # supporting reads already counted, so now we can lump in the breakpoint-mapping reads as well
    trna.ind <- which( trna.ids %in% unique( c( ids.support, ids.maybe ) ) )
    # identify nwes reads that should be remapped
    ids.maybe <- suppressWarnings( system2( command='samtools',
                                   args=c( 'view',
                                           '-@ 12',
                                           'tmp/nwes/recal.bam',
                                           paste(paste(gene1.chr,bpt1.L,sep=':'),bpt1.R,sep='-'),
                                           paste(paste(gene2.chr,bpt2.L,sep=':'),bpt2.R,sep='-'),
                                           '| cut -f1 | sort -u' ),
                                   stdout=T, stderr=T ) )
    nwes.ind <- which( nwes.ids %in% unique( ids.maybe ) )
    # identify twes reads that should be remapped
    ids.maybe <- suppressWarnings( system2( command='samtools',
                                   args=c( 'view',
                                           '-@ 12',
                                           'tmp/twes/recal.bam',
                                           paste(paste(gene1.chr,bpt1.L,sep=':'),bpt1.R,sep='-'),
                                           paste(paste(gene2.chr,bpt2.L,sep=':'),bpt2.R,sep='-'),
                                           '| cut -f1 | sort -u' ),
                                   stdout=T, stderr=T ) )
    twes.ind <- which( twes.ids %in% unique( ids.maybe ) )
    # ... base directory for output files, number of supporting read pairs prepended
    bn <- paste( 'tmp/trna/fusionqc/numsrp',
                 paste( sprintf( "%04d", results$support_count[i] ),
                        pair.underscores,
                        paste( gene1.chr, gene1.bpt, gene1.std, sep='_' ),
                        paste( gene2.chr, gene2.bpt, gene2.std, sep='_' ),
                        sep='__' ),
                 sep='' )
    # write test reference and annotation tracks
    dir.create( path=bn )
    # ... ref fasta, indexed
    testref.seq <- c( gene1.seq, gene2.seq, fusion.seq )
    writeXStringSet( testref.seq, paste( bn, '/testref.fa', sep='' ) )
    indexFa( paste( bn, '/testref.fa', sep='' ) )  # create .fai
    buildindex( basename=paste( bn, '/testref', sep='' ), reference=paste( bn, '/testref.fa', sep='' ) )  # index for subjunc
    # ... annotation of fusion parts
    fuseparts.ann <- suppressWarnings( c( fuse1.gr, fuse2.gr, fused.gr ) )
    export( fuseparts.ann, paste( bn, '/fusion_parts.gtf', sep='' ) )
    # ... annotation of gene parts
    testref.ann <- GRanges( seqnames=c( seqnames( gene1.ann ), seqnames( gene2.ann ), seqnames( fusion.ann ) ),
                            strand=c( strand( gene1.ann ), strand( gene2.ann ), strand( fusion.ann ) ),
                            ranges=c( ranges( gene1.ann ), ranges( gene2.ann ), ranges( fusion.ann ) ) )
    mcols( testref.ann ) <- rbind( mcols( gene1.ann ), mcols( gene2.ann ), mcols( fusion.ann ) )
    export( testref.ann, paste( bn, '/gene_parts.gtf', sep='' ) )
    # pull and realign reads, create BAMs
    # ... should we consider {subjunc} detectSV=T and reportAllJunctions=T for fusion / altsplice calling??
    # ... start track_config.json text block:
    config <- list( list( name = 'fusion_parts',
                          url = 'fusion_parts.gtf',
                          height = 75,
                          colorBy = 'gene_name',
                          colorTable = list( fusion_partner_1 = 'darkgreen',
                                             fusion_partner_2 = 'blueviolet' ) ),
                    list( name = 'gene_parts',
                          url = 'gene_parts.gtf',
                          height = 250 ) )
    # ... trna:
    if ( length( trna.ind ) > 0 ) {
      temp.R1 <- trna.R1[ trna.ind ]
      temp.R2 <- trna.R2[ trna.ind ]
      writeFastq( temp.R1, file=paste( bn, '/trna.support_R1.fq.gz', sep='' ) )
      writeFastq( temp.R2, file=paste( bn, '/trna.support_R2.fq.gz', sep='' ) )
      # ... ... align
      subjunc( index=paste( bn, '/testref', sep='' ),
               readfile1=paste( bn, '/trna.support_R1.fq.gz', sep='' ),
               readfile2=paste( bn, '/trna.support_R2.fq.gz', sep='' ),
               sortReadsByCoordinates=T, useAnnotation=T, isGTF=T,
               annot.ext=paste( bn, '/gene_parts.gtf', sep='' ),
               minFragLength=35, maxFragLength=1500,
               nthreads=2,
               output_file=paste( bn, '/trna.support.bam', sep='' ) )
      # ... ... create MQ > 0 bam
      temp <- suppressWarnings( system2( command='samtools',
                                         args=c( 'view',
                                                 '-q 1',
                                                 '-b',
                                                 paste('-o ',bn,'/trna.support.MQgt0.bam',sep=''),
                                                 paste(bn,'/trna.support.bam',sep='') ),
                                         stdout=T, stderr=T ) )
      # ... ... and index it
      temp <- suppressWarnings( system2( command='samtools',
                                         args=c( 'index',
                                                 paste(bn,'/trna.support.MQgt0.bam',sep='') ),
                                         stdout=T, stderr=T ) )
      # ... ... track config
      config.trna <- list( list( name = 'trna.support.MQgt0',
                                 url = 'trna.support.MQgt0.bam',
                                 viewAsPairs = 'true' ),
                           list( name = 'trna.support',
                                 url = 'trna.support.bam',
                                 viewAsPairs = 'true' ) )
      config <- c( config, config.trna )
    }
    # ... twes:
    if ( length( twes.ind ) > 0 ) {
      temp.R1 <- twes.R1[ twes.ind ]
      temp.R2 <- twes.R2[ twes.ind ]
      writeFastq( temp.R1, file=paste( bn, '/twes.support_R1.fq.gz', sep='' ) )
      writeFastq( temp.R2, file=paste( bn, '/twes.support_R2.fq.gz', sep='' ) )
      # ... ... align
      align( index=paste( bn, '/testref', sep='' ),
             readfile1=paste( bn, '/twes.support_R1.fq.gz', sep='' ),
             readfile2=paste( bn, '/twes.support_R2.fq.gz', sep='' ),
             type='dna', sortReadsByCoordinates=T,
             minFragLength=35, maxFragLength=1500,
             nthreads=4,
             output_file=paste( bn, '/twes.support.bam', sep='' ) )
      # ... ... create MQ > 0 bam
      temp <- suppressWarnings( system2( command='samtools',
                                         args=c( 'view',
                                                 '-q 1',
                                                 '-b',
                                                 paste('-o ',bn,'/twes.support.MQgt0.bam',sep=''),
                                                 paste(bn,'/twes.support.bam',sep='') ),
                                         stdout=T, stderr=T ) )
      # ... ... and index it
      temp <- suppressWarnings( system2( command='samtools',
                                         args=c( 'index',
                                                 paste(bn,'/twes.support.MQgt0.bam',sep='') ),
                                         stdout=T, stderr=T ) )
      # ... ... track config
      config.twes <- list( list( name = 'twes.support.MQgt0',
                                 url = 'twes.support.MQgt0.bam',
                                 viewAsPairs = 'true' ),
                           list( name = 'twes.support',
                                 url = 'twes.support.bam',
                                 viewAsPairs = 'true' ) )
      config <- c( config, config.twes )
    }
    # ... nwes:
    if ( length( nwes.ind ) > 0 ) {
      temp.R1 <- nwes.R1[ nwes.ind ]
      temp.R2 <- nwes.R2[ nwes.ind ]
      writeFastq( temp.R1, file=paste( bn, '/nwes.support_R1.fq.gz', sep='' ) )
      writeFastq( temp.R2, file=paste( bn, '/nwes.support_R2.fq.gz', sep='' ) )
      # ... ... align
      align( index=paste( bn, '/testref', sep='' ),
             readfile1=paste( bn, '/nwes.support_R1.fq.gz', sep='' ),
             readfile2=paste( bn, '/nwes.support_R2.fq.gz', sep='' ),
             type='dna', sortReadsByCoordinates=T,
             minFragLength=35, maxFragLength=1500,
             nthreads=2,
             output_file=paste( bn, '/nwes.support.bam', sep='' ) )
      # ... ... create MQ > 0 bam
      temp <- suppressWarnings( system2( command='samtools',
                                         args=c( 'view',
                                                 '-q 1',
                                                 '-b',
                                                 paste('-o ',bn,'/nwes.support.MQgt0.bam',sep=''),
                                                 paste(bn,'/nwes.support.bam',sep='') ),
                                         stdout=T, stderr=T ) )
      # ... ... and index it
      temp <- suppressWarnings( system2( command='samtools',
                                         args=c( 'index',
                                                 paste(bn,'/nwes.support.MQgt0.bam',sep='') ),
                                         stdout=T, stderr=T ) )
      # ... ... track config
      config.nwes <- list( list( name = 'nwes.support.MQgt0',
                                 url = 'nwes.support.MQgt0.bam',
                                 viewAsPairs = 'true' ),
                           list( name = 'nwes.support',
                                 url = 'nwes.support.bam',
                                 viewAsPairs = 'true' ) )
      config <- c( config, config.nwes )
    }
    # ... write track config:
    write_json( x = config, auto_unbox = T,
                path = paste(bn,'/track_config.json',sep='') )
    # delete large shortread index files after all remapping completed:
    unlink( paste( bn, '/testref.00.*', sep='' ) )
    # create html IGV Report
    # ... sites file
    temp <- data.frame( sequence=c( paste('gene1',gene1.name,sep='_'),
                                    paste('gene2',gene2.name,sep='_'),
                                    pair.underscores ),
                        begin=c( 1, 1, 1 ),
                        end=c( nchar( gene1.seq ), nchar( gene2.seq ), nchar( fusion.seq ) ) )
    write.table( temp, file=paste(bn,'/sites.tsv',sep=''), col.names=T, row.names=F, quote=F, sep="\t" )
    # ... create the html report
    # (add '--translate-sequence-track' once new tagged release is avbl via bioconda)
    #file.copy( '/opt/track_config.json', paste(bn,'/',sep='') )
    temp <- suppressWarnings( system2( command='create_report',
                                       args=c( paste(bn,'/sites.tsv',sep=''),
                                               paste(bn,'/testref.fa',sep=''),
                                               '--begin 2',
                                               '--end 3',
                                               '--sequence 1',
                                               paste('--track-config',paste(bn,'/track_config.json',sep=''),sep=' '),
                                               paste('--output',paste(bn,'/igvreport.html',sep=''),sep=' '),
                                               '--standalone'
                                             ), stdout=T, stderr=T ) )
    # overwrite results with additional entry (potentially allows for restarting without progress loss, if error)
    write.table( results[,-11], file='tmp/trna/fusionqc/results.tsv', quote=F, row.names=F, sep="\t" )
  }
}


