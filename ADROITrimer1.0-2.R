#==================================================
#==================================================
#==================================================
# Libraries

library( bio3d )



#==================================================
#==================================================
#==================================================
# Functions

#==================================================
# Check if x is a letter
is.letter <- function( x )
{
  grepl( "[[:alpha:]]",
         x )
}

#==================================================
# Main ADROITrimer function
ADROITrimer <- function( file.reference.alignment,
                         file.repair.sequence,
                         subtype,
                         exe.muscle,
                         input_prefix = paste( "input",
                                               round( as.numeric( Sys.time() ,digits = 10 ) ),
                                               sep = "_" ),
                         output.prefix = paste( "output",
                                               round( as.numeric( Sys.time() ,digits = 10 ) ),
                                               sep = "_" ) )
{
  #==================================================
  # Set variables
  print( "Set variables" )
  
  file.reference.alignment.subtype <- paste( unlist( strsplit( file.reference.alignment,
                                                               "\\.fasta" ) ),
                                             "_subtype",
                                             subtype,
                                             ".fasta",
                                             sep = "" )
  file.reference.alignment.subtype.repair <- paste( unlist( strsplit( file.reference.alignment.subtype,
                                                                      "\\.fasta" ) ),
                                                    "_",
                                                    input_prefix,
                                                    ".fasta",
                                                    sep = "" )
  file.reference.alignment.subtype.repair.muscle <- paste( unlist( strsplit( file.reference.alignment.subtype,
                                                                             "\\.fasta" ) ),
                                                           "_",
                                                           input_prefix,
                                                           "_muscle.fasta",
                                                           sep = "" )
  file.reference.alignment.subtype.repair.muscle.no.hxb2 <- paste( unlist( strsplit( file.reference.alignment.subtype,
                                                                                     "\\.fasta" ) ),
                                                                   "_",
                                                                   input_prefix,
                                                                   "_muscle_noHXB2.fasta",
                                                                   sep = "" )
  file.new.aln <- paste( unlist( strsplit( file.reference.alignment.subtype,
                                           "\\.fasta" ) ),
                         "_",
                         output.prefix,
                         ".fasta",
                         sep = "" )
  
  #==================================================
  # Load raw HIV alignment and select for subtype of interest
  print( "Load raw HIV alignment and select for subtype of interest" )
  
  aln.ref <- read.fasta( file.reference.alignment )
  aln.ref.ali <- aln.ref$ali
  aln.ref.id <- aln.ref$id
  
  vec.subtypes <- unlist( lapply( aln.ref.id,
                                  function( id )
                                  {
                                    return( unlist( strsplit( id,
                                                              "\\." ) )[ 1 ] )
                                  } ) )
  
  # Clean subtypes
  vec.subtypes.clean <- NULL
  for( type in vec.subtypes )
  {
    vec.char.UPPER <- NULL
    for( char in unlist( strsplit( type,
                                   "" ) ) )
    {
      if( char == toupper( char ) &
          is.letter( char ) )
      {
        vec.char.UPPER <- c( vec.char.UPPER,
                             char )
      }
    }
    vec.subtypes.clean <- c( vec.subtypes.clean,
                             paste( unique( vec.char.UPPER ),
                                    collapse = "" ) )
  }
  
  #==================================================
  # Identify if sufficient number (>=50) of subtypes available
  print( "Identify if sufficient number (>=50) of subtypes available" )
  
  vec.subtypes.ind <- which( vec.subtypes.clean == subtype )
  if( length( vec.subtypes.ind ) >= 50 )
  {
    # Does this subtype-only alignment exist?
    if( !( file.exists( file.reference.alignment.subtype ) ) )
    {
      write.fasta( ids = aln.ref.id[ c( 1, vec.subtypes.ind ) ],
                   seqs = aln.ref.ali[ c( 1, vec.subtypes.ind ), ],
                   file = file.reference.alignment.subtype )
    }
  } else
  {
    file.reference.alignment.subtype <- file.reference.alignment
  }
  
  # Load subtype only sequences
  aln.subtype <- read.fasta( file.reference.alignment.subtype )
  aln.subtype.ali <- aln.subtype$ali
  aln.subtype.id <- aln.subtype$id
  
  #==================================================
  # Load test sequence
  print( "Load test sequence" )
  
  aln.test <- read.fasta( file.repair.sequence )
  aln.test.ali <- aln.test$ali
  aln.test.id <- aln.test$id
  
  #==================================================
  # Align test sequence with subtype alignment
  print( "Align test sequence with subtype alignment" )
  
  aln.subtype.test.ali <- rbind( aln.subtype.ali,
                                 c( aln.test.ali,
                                    rep( "-",
                                         ncol( aln.subtype.ali ) - length( aln.test.ali ) ) ) )
  aln.subtype.test.id <- c( aln.subtype.id,
                            aln.test.id )
  write.fasta( ids = aln.subtype.test.id,
               seqs = aln.subtype.test.ali,
               file = file.reference.alignment.subtype.repair )
  
  # Align using muscle with 2 iterations
  system( paste( exe.muscle,
                 "-in",
                 file.reference.alignment.subtype.repair,
                 "-out",
                 file.reference.alignment.subtype.repair.muscle,
                 "-maxiters",
                 "2" ) )
  
  #==================================================
  # Load new subtype alignment including HXB2 and test sequence
  print( "Load new subtype alignment including HXB2 and test sequence" )
  
  aln.subtype.test.muscle <- read.fasta( file.reference.alignment.subtype.repair.muscle )
  aln.subtype.test.muscle.ali <- aln.subtype.test.muscle$ali
  aln.subtype.test.muscle.id <- aln.subtype.test.muscle$id
  
  
  # Remove (initial) HXB2 sequences from alignment and save separately
  var.ind.hxb2 <- which( grepl( "HXB2",
                                aln.subtype.test.muscle.id ) )
  vec.seq.hxb2 <- aln.subtype.test.muscle.ali[ var.ind.hxb2[ 1 ], ]
  
  aln.subtype.test.muscle.ali.no.hxb2 <- aln.subtype.test.muscle.ali[ -var.ind.hxb2[ 1 ], ]
  aln.subtype.test.muscle.id.no.hxb2 <- aln.subtype.test.muscle.id[ -var.ind.hxb2[ 1 ] ]
  
  # Save aligned subtype alignment without reference sequence HXB2
  write.fasta( ids = aln.subtype.test.muscle.id.no.hxb2,
               seqs = aln.subtype.test.muscle.ali.no.hxb2,
               file = file.reference.alignment.subtype.repair.muscle.no.hxb2 )
  
  aln.subtype.test.muscle.no.hxb2 <- read.fasta( file.reference.alignment.subtype.repair.muscle.no.hxb2 )
  aln.subtype.test.muscle.no.hxb2.ali <- aln.subtype.test.muscle.no.hxb2$ali
  aln.subtype.test.muscle.no.hxb2.id <- aln.subtype.test.muscle.no.hxb2$id
  
  # Extract aligned test sequence
  var.ind.test <- which( aln.test.id == aln.subtype.test.muscle.id )
  vec.seq.test.aligned <- aln.subtype.test.muscle.ali[ var.ind.test, ]
  
  
  
  #==================================================
  # Generate initial data frame
  print( "Generate initial data frame" )
  
  vec.inserts <- c( LETTERS,
                    letters )
  df.no_aln.no_hxb2.insert_hxb2.hxb2_aa <- NULL
  ind.hxb2 <- 0
  ind.letters <- 0
  for( i in 1:ncol( aln.subtype.test.muscle.no.hxb2.ali ) )
  {
    if( vec.seq.hxb2[ i ] != "-" )
    {
      ind.hxb2 <- ind.hxb2 + 1
      insert <- " "
      ind.letters <- 0
      # no.hxb2 <- ind.hxb2
    } else
    {
      ind.letters <- ind.letters + 1
      # no.hxb2 <- paste( ind.hxb2,
      #                   vec.inserts[ ind.letters ],
      #                   sep = "" )
      insert <- vec.inserts[ ind.letters ]
    }
    df.no_aln.no_hxb2.insert_hxb2.hxb2_aa <- rbind( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa,
                                                    c( i,
                                                       ind.hxb2,
                                                       insert,
                                                       vec.seq.hxb2[ i ],
                                                       vec.seq.test.aligned[ i ] ) )
  }
  df.no_aln.no_hxb2.insert_hxb2.hxb2_aa <- data.frame( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa )
  colnames( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa ) <- c( "no_aln",
                                                          "no_hxb2",
                                                          "insert_hxb2",
                                                          "aa_hxb2",
                                                          "aa_newSequence" )
  
  #==================================================
  # Calculate consensus sequence with default (0.6) cutoff
  print( "Calculate consensus sequence with default (0.6) cutoff" )
  
  con <- consensus( aln.subtype.test.muscle.no.hxb2 )
  df.AA.freq <- data.frame( t( con$freq ) )
  colnames( df.AA.freq ) <- colnames( t( con$freq ) )
  vec.con.names <- colnames( t( con$freq ) )
  
  #==================================================
  # Data frame containing consensus frequencies
  print( "Data frame containing consensus frequencies" )
  
  aa_con <- con$seq
  df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq <- cbind( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa,
                                                          aa_con,
                                                          df.AA.freq )
  
  #==================================================
  # Upload PDB
  pdb <- read.pdb( "5fyl_trimer.pdb" )
  pdb.atom <- pdb$atom
  pdb.atom.protomer <- pdb.atom[ pdb.atom$chain %in% c( "A",
                                                        "B" ) &
                                   pdb.atom$elety == "CA", ]
  
  # Upload RSAs
  lines.rsa <- readLines( "5fyl_trimer.rsa" )
  ind.res.hem <- c( which( grepl("^RES ",
                                 lines.rsa ) ),
                    which( grepl("^HEM ",
                                 lines.rsa ) ) )
  lines.rsa.trunc <- lines.rsa[ ind.res.hem ]
  
  df.rsa <- NULL
  for( line.rsa in lines.rsa.trunc )
  {
    line <- unlist( strsplit( line.rsa,
                              "" ) )
    df.rsa <- rbind( df.rsa,
                     c( paste( line[ 1:3 ],
                               collapse = "" ),
                        paste( line[ 5:7 ],
                               collapse = "" ),
                        line[ 9 ],
                        as.numeric( paste( line[ 10:13 ],
                                           collapse = "" ) ),
                        as.numeric( paste( line[ 37:41 ],
                                           collapse = "" ) ) ) )
  }
  
  df.rsa <- data.frame( df.rsa )
  colnames( df.rsa ) <- c( "atom",
                           "aa",
                           "chain",
                           "resno",
                           "sasa" )
  df.rsa.protomer <- df.rsa[ df.rsa$chain %in% c( "A",
                                                  "B" ), ]
  df.rsa.protomer$sasa <- as.numeric( as.vector( df.rsa.protomer$sasa ) )
  
  #==================================================
  # Merge RSA with data frame
  df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa <- cbind( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq,
                                                              rep( NA,
                                                                   nrow( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq ) ) )
  colnames( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa ) <- c( colnames( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq ),
                                                                      "RSA" )
  
  df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$no_hxb2 <- as.numeric( as.vector( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$no_hxb2 ) )
  df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$insert_hxb2 <- as.vector( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$insert_hxb2 )
  df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$aa_hxb2 <- as.vector( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$aa_hxb2 )
  
  for( i in 1:nrow( pdb.atom.protomer ) )
  {
    #==================================================
    # Manually ignore positions 321 and 322
    if( pdb.atom.protomer[ i, ]$resno %in% 321:322 )
    {
      next
    }
    #==================================================
    if( is.na( pdb.atom.protomer[ i, ]$insert ) )
    {
      df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa[ df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$no_hxb2 == pdb.atom.protomer[ i, ]$resno &
                                                           df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$insert_hxb2 == " " &
                                                           !( is.na( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$insert_hxb2 ) ), ]$RSA <- df.rsa.protomer[ i, ]$sasa
    } else
    {
      df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa[ df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$no_hxb2 == pdb.atom.protomer[ i, ]$resno &
                                                           df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$insert_hxb2 == pdb.atom.protomer[ i, ]$insert, ]$RSA <- df.rsa.protomer[ i, ]$sasa
    }
  }
  df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa <- data.frame( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa )
  df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$aa_newSequence <- as.vector( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$aa_newSequence )
  df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$aa_con <- as.vector( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$aa_con )
  df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$no_hxb2 <- as.numeric( as.vector( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$no_hxb2 ) )
  
  #==================================================
  # Repair sequences
  # 1. Mutate residue to consensus residue if the prevalence of WT residue is less than 2%
  # 2. Mutate residue to consensus residue if the prevalence of WT residue is between 2% to 75% AND the buried surface area is less than 30%?
  print( "Repair sequence" )
  
  vec.seq.repair <- NULL
  # 10% of RSA
  # var.rsa.10percent <- 0.1 * max( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$RSA,
  #                                 na.rm = TRUE )
  # 30% of RSA
  var.rsa.30percent <- 0.3 * max( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa$RSA,
                                  na.rm = TRUE )
  vec.aa.new.seq.freq <- rep( NA,
                              nrow( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa ) )
  
  for( i in 1:nrow( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa ) )
  {
    aa.new.seq <- df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa[ i, ]$aa_newSequence
    
    #==================================================
    # Manually ignore positions 321 and 322
    if( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa[ i, ]$no_hxb2 %in% 321:322 )
    {
      vec.seq.repair <- c( vec.seq.repair,
                           aa.new.seq )
      next
    }
    #==================================================
    if( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa[ i, ]$no_hxb2 %in% 508:511 )
    {
      vec.seq.repair <- c( vec.seq.repair,
                           aa.new.seq )
    } else if( aa.new.seq == "-" )
    {
      vec.seq.repair <- c( vec.seq.repair,
                           aa.new.seq )
      
    } else
    {
      j <- which( vec.con.names == aa.new.seq )
      vec.aa.new.seq.freq[ i ] <- df.AA.freq[ i, j ]
      if( df.AA.freq[ i, j ] < 0.02 )
      {
        if( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa[ i, ]$aa_con == "-" )
        {
          vec.seq.repair <- c( vec.seq.repair,
                               aa.new.seq )
        } else
        {
          vec.seq.repair <- c( vec.seq.repair,
                               df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa[ i, ]$aa_con ) 
        }
      } else if( df.AA.freq[ i, j ] >= 0.02 &
                 df.AA.freq[ i, j ] < 0.075 )
      {
        if( is.na( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa[ i, ]$RSA ) )
        {
          vec.seq.repair <- c( vec.seq.repair,
                               aa.new.seq )
        } else
        {
          if( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa[ i, ]$RSA < var.rsa.30percent )
          {
            if( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa[ i, ]$aa_con == "-" )
            {
              vec.seq.repair <- c( vec.seq.repair,
                                   aa.new.seq )
            } else
            {
              vec.seq.repair <- c( vec.seq.repair,
                                   df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa[ i, ]$aa_con ) 
            }
          } else
          {
            vec.seq.repair <- c( vec.seq.repair,
                                 aa.new.seq )
          }
        }
      } else
      {
        vec.seq.repair <- c( vec.seq.repair,
                             aa.new.seq )
      }
    }
  }
  
  
  df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq <- cbind( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa,
                                                                          vec.aa.new.seq.freq )
  colnames( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq ) <- c( colnames( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa ),
                                                                                  "aa_new_seq_freq" )
  #==================================================
  # Data frame containing raw sequence and repaired version
  df.newseq.repairseq <- data.frame( cbind( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq$aa_newSequence,
                                            vec.seq.repair ) )
  colnames( df.newseq.repairseq ) <- c( "new_seq",
                                        "repair_seq" )
  df.newseq.repairseq$new_seq <- as.vector( df.newseq.repairseq$new_seq )
  df.newseq.repairseq$repair_seq <- as.vector( df.newseq.repairseq$repair_seq )
  
  
  #==================================================
  # Add two additional columns with repaired sequence and mutation performed
  df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation <- data.frame( cbind( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq,
                                                                                                           vec.seq.repair,
                                                                                                           rep( NA,
                                                                                                                nrow( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq ) ) ) )
  colnames( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation ) <- c( colnames( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq ),
                                                                                                       "repaired_aa",
                                                                                                       "mutation" )
  df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation$repaired_aa <- as.vector( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation$repaired_aa )
  df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation$no_hxb2 <- as.vector( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation$no_hxb2 )
  df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation$insert_hxb2 <- as.vector( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation$insert_hxb2 )
  
  vec.mutation.list <- NULL
  for( i in 1:nrow( df.newseq.repairseq ) )
  {
    if( df.newseq.repairseq[ i, 1 ] != df.newseq.repairseq[ i , 2 ] )
    {
      if( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation[ i, ]$insert_hxb2 == " " )
      {
        df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation[ i, ]$mutation <- paste( df.newseq.repairseq[ i, 1 ],
                                                                                                                    df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation[ i, ]$no_hxb2,
                                                                                                                    df.newseq.repairseq[ i, 2 ],
                                                                                                                    sep = "" )
        vec.mutation.list <- c( vec.mutation.list,
                                paste( df.newseq.repairseq[ i, 1 ],
                                       df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation[ i, ]$no_hxb2,
                                       df.newseq.repairseq[ i, 2 ],
                                       sep = "" ) )
      } else
      {
        df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation[ i, ]$mutation <- paste( df.newseq.repairseq[ i, 1 ],
                                                                                                                    df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation[ i, ]$no_hxb2,
                                                                                                                    df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation[ i, ]$insert_hxb2,
                                                                                                                    df.newseq.repairseq[ i, 2 ],
                                                                                                                    sep = "" )
        vec.mutation.list <- c( vec.mutation.list,
                                paste( df.newseq.repairseq[ i, 1 ],
                                       df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation[ i, ]$no_hxb2,
                                       df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation[ i, ]$insert_hxb2,
                                       df.newseq.repairseq[ i, 2 ],
                                       sep = "" ) )
      }
    }
  }
  
  
  
  #==================================================
  # Save final CSV file
  write.csv( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation,
             file = paste( output.prefix,
                           "_FullAnalysis.csv",
                           sep = "" ),
             row.names = FALSE,
             quote = FALSE )
  
  # Save file with just the mutations
  write( vec.mutation.list,
         file = paste( output.prefix,
                       "_repaired_mutList.txt",
                       sep = "" ),
         ncolumns = 1 )
  
  #==================================================
  # Save new repaired sequences as FASTA file
  vec.seq.repair.no.gaps <- vec.seq.repair[ vec.seq.repair != "-" ]
  write.fasta( seqs = vec.seq.repair.no.gaps,
               id = output.prefix,
               file = paste( output.prefix,
                             "_repaired.fasta",
                             sep = "" ) )
  
  #==================================================
  #==================================================
  # Additional mutations
  print( "Additional mutations" )
  
  df.no_aln.no_hxb2.insert_hxb2.repaired_aa <- df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation[ , c( 1:3, 
                                                                                                                                        ncol( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation ) - 1,
                                                                                                                                        ncol( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation ) ) ]
  # (i)   Add SOSIP mutation
  # (ii)  Add DS
  # (iii) Add 6R before gp41
  # (iv)  Add Janssen mutations
  # (v)   Truncate after 664
  # (vi)  Add Nt and Ct termini
  vec.final.sequence <- NULL
  df.tmp <- NULL
  for( i in 1:nrow( df.no_aln.no_hxb2.insert_hxb2.repaired_aa ) )
  {
    # Remove N-terminal part
    if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] %in% 1:34 )
    {
      next
    } # Truncate after 664
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 665 )
    {
      break
    } # Add 6R
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 508 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "R" )
      df.tmp <- rbind( df.tmp,
                       c( 508, "6R" ) )
    } # Add 6R
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 509 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "R" )
      df.tmp <- rbind( df.tmp,
                       c( 509, "6R" ) )
    } # Add 6R
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 510 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "R" )
      df.tmp <- rbind( df.tmp,
                       c( 510, "6R" ) )
    } # Add 6R
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 511 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "RRR" )
      df.tmp <- rbind( df.tmp,
                       c( 511, "6R" ),
                       c( "511a", "6R" ),
                       c( "511b", "6R" ) )
    } # Add DS
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] %in% c( 201,
                                                                            433 ) &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "C" )
      df.tmp <- rbind( df.tmp,
                       c( 201, "DS" ),
                       c( 433, "DS" ) )
    } # Add SOSIP
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] %in% c( 501,
                                                                            605 ) &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "C" )
      df.tmp <- rbind( df.tmp,
                       c( 501, "SOSIP" ),
                       c( 605, "SOSIP" ) )
    } # Add SOSIP
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 559 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "P" )
      df.tmp <- rbind( df.tmp,
                       c( 559, "SOSIP" ) )
    } # Add Janssen mutation
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 535 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "N" )
      df.tmp <- rbind( df.tmp,
                       c( 535, "Stabilization" ) )
    } # Add Janssen mutation
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 556 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "P" )
      df.tmp <- rbind( df.tmp,
                       c( 556, "Stabilization" ) )
    } # Add Janssen mutation
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 588 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "E" )
      df.tmp <- rbind( df.tmp,
                       c( 588, "Stabilization" ) )
    } # Add Janssen mutation
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 589 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "V" )
      df.tmp <- rbind( df.tmp,
                       c( 589, "Stabilization" ) )
    } # Add Janssen mutation
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 651 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "F" )
      df.tmp <- rbind( df.tmp,
                       c( 651, "Stabilization" ) )
    } # Add Janssen mutation
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 655 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "I" )
      df.tmp <- rbind( df.tmp,
                       c( 655, "Stabilization" ) )
    } # Add Janssen mutation
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 658 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "V" )
      df.tmp <- rbind( df.tmp,
                       c( 658, "Stabilization" ) )
    } 
    # # Add 368R mutation (Mangai's request)
    # else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 368 &
    #          df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    # {
    #   vec.final.sequence <- c( vec.final.sequence,
    #                            "R" )
    # }
    else
    {
      vec.final.sequence <- c( vec.final.sequence,
                               df.no_aln.no_hxb2.insert_hxb2.repaired_aa$repaired_aa[ i ] )

      if( !( is.na( df.no_aln.no_hxb2.insert_hxb2.repaired_aa[ i, 5 ] ) ) )
      {
        if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa[ i, 3 ] == " " )
        {
          df.tmp <- rbind( df.tmp,
                           c( df.no_aln.no_hxb2.insert_hxb2.repaired_aa[ i, 2 ],
                              paste( "Repair",
                                     df.no_aln.no_hxb2.insert_hxb2.repaired_aa[ i, 5 ],
                                     sep = "_" ) ) )
        } else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa[ i, 3 ] != " " )
        {
          df.tmp <- rbind( df.tmp,
                           c( paste( df.no_aln.no_hxb2.insert_hxb2.repaired_aa[ i, 2 ],
                                     df.no_aln.no_hxb2.insert_hxb2.repaired_aa[ i, 3 ],
                                     sep = "_" ),
                              paste( "Repair",
                                     df.no_aln.no_hxb2.insert_hxb2.repaired_aa[ i, 5 ],
                                     sep = "_" ) ) )
        }
      }
    }
  }
  
  # Add N- and C-terminal tags
  vec.final.sequence.Nt.Ct <- c( unlist( strsplit( "MPMGSLQPLATLYLLGMLVASVLAQQAENL",
                                                   "" ) ),
                                 vec.final.sequence,
                                 unlist( strsplit( "GSAPTKAKRRVVQREKR",
                                                   "" ) ) )
  # vec.final.sequence.Nt.Ct <- c( vec.final.sequence )
  
  #==================================================
  # Save repaired sequence with additional mutations
  vec.final.sequence.Nt.Ct <- vec.final.sequence.Nt.Ct[ vec.final.sequence.Nt.Ct != "-" ]
  write.fasta( seqs = vec.final.sequence.Nt.Ct,
               id = output.prefix,
               file = paste( output.prefix,
                             "_DS-SOSIP-RnS_NtCt-tags.fasta",
                             sep = "" ) )
  
  df.tmp <- data.frame( df.tmp )
  colnames( df.tmp ) <- c( "HXB2",
                           "Mutation" )
  write.csv( df.tmp,
             file = "Mutations.csv",
             row.names = FALSE,
             quote = FALSE )
  
  
  #==================================================
  #==================================================
  # Additional mutations
  print( "Additional mutations" )
  
  df.no_aln.no_hxb2.insert_hxb2.repaired_aa <- df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation[ , c( 1:3, 
                                                                                                                                        ncol( df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation ) - 1 ) ]
  # (i)   Add SOSIP mutation
  # (ii)  Add DS
  # (iii) Add 6R before gp41
  # (iv)  Add Janssen mutations
  # (v)   Truncate after 664
  vec.final.sequence <- NULL
  for( i in 1:nrow( df.no_aln.no_hxb2.insert_hxb2.repaired_aa ) )
  {
    # Remove N-terminal part
    if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] %in% 1:34 )
    {
      next
    } # Truncate after 664
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 665 )
    {
      break
    } # Add 6R
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 508 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "R" )
    } # Add 6R
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 509 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "R" )
    } # Add 6R
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 510 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "R" )
    } # Add 6R
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 511 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "RRR" )
    } # Add DS
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] %in% c( 201,
                                                                            433 ) &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "C" )
    } # Add SOSIP
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] %in% c( 501,
                                                                            605 ) &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "C" )
    } # Add SOSIP
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 559 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "P" )
    } # Add Janssen mutation
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 535 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "N" )
    } # Add Janssen mutation
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 556 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "P" )
    } # Add Janssen mutation
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 588 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "E" )
    } # Add Janssen mutation
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 589 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "V" )
    } # Add Janssen mutation
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 651 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "F" )
    } # Add Janssen mutation
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 655 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "I" )
    } # Add Janssen mutation
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 658 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "V" )
    } 
    # # Add 368R mutation (Mangai's request)
    # else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 368 &
    #          df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    # {
    #   vec.final.sequence <- c( vec.final.sequence,
    #                            "R" )
    # }
    else
    {
      vec.final.sequence <- c( vec.final.sequence,
                               df.no_aln.no_hxb2.insert_hxb2.repaired_aa$repaired_aa[ i ] )
    }
  }
  
  # Add N- and C-terminal tags
  # vec.final.sequence.Nt.Ct <- c( unlist( strsplit( "MPMGSLQPLATLYLLGMLVASVLAQQAENL",
  #                                                  "" ) ),
  #                                vec.final.sequence,
  #                                unlist( strsplit( "GSAPTKAKRRVVQREKR",
  #                                                  "" ) ) )
  vec.final.sequence.Nt.Ct <- c( unlist( strsplit( "AENL",
                                                   "" ) ),
                                 vec.final.sequence )
  
  #==================================================
  # Save repaired sequence with additional mutations
  vec.final.sequence.Nt.Ct <- vec.final.sequence.Nt.Ct[ vec.final.sequence.Nt.Ct != "-" ]
  write.fasta( seqs = vec.final.sequence.Nt.Ct,
               id = output.prefix,
               file = paste( output.prefix,
                             "_DS-SOSIP-RnS.fasta",
                             sep = "" ) )
  
  #==================================================
  #==================================================
  # Save DS-SOSIP only
  
  df.no_aln.no_hxb2.insert_hxb2.repaired_aa <- df.no_aln.no_hxb2.insert_hxb2.hxb2_aa.AA_freq.rsa.aa_new_freq.repaired_aa.mutation[ , c( 1:3, 5 ) ]
  # (i)   Add SOSIP mutation
  # (ii)  Add DS
  # (iii) Add 6R before gp41
  # (iv)   Truncate after 664
  # (v)  Add Nt and Ct termini
  vec.final.sequence <- NULL
  for( i in 1:nrow( df.no_aln.no_hxb2.insert_hxb2.repaired_aa ) )
  {
    # Remove N-terminal part
    if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] %in% 1:34 )
    {
      next
    } # Truncate after 664
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 665 )
    {
      break
    } # Add 6R
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 508 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "R" )
    } # Add 6R
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 509 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "R" )
    } # Add 6R
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 510 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "R" )
    } # Add 6R
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 511 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "RRR" )
    } # Add DS
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] %in% c( 201,
                                                                            433 ) &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "C" )
    } # Add SOSIP
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] %in% c( 501,
                                                                            605 ) &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "C" )
    } # Add SOSIP
    else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 559 &
             df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    {
      vec.final.sequence <- c( vec.final.sequence,
                               "P" )
    } # Add Janssen mutation
    # else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 535 &
    #          df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    # {
    #   vec.final.sequence <- c( vec.final.sequence,
    #                            "N" )
    # } # Add Janssen mutation
    # else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 556 &
    #          df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    # {
    #   vec.final.sequence <- c( vec.final.sequence,
    #                            "P" )
    # } # Add Janssen mutation
    # else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 588 &
    #          df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    # {
    #   vec.final.sequence <- c( vec.final.sequence,
    #                            "E" )
    # } # Add Janssen mutation
    # else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 589 &
    #          df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    # {
    #   vec.final.sequence <- c( vec.final.sequence,
    #                            "V" )
    # } # Add Janssen mutation
    # else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 651 &
    #          df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    # {
    #   vec.final.sequence <- c( vec.final.sequence,
    #                            "F" )
    # } # Add Janssen mutation
    # else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 655 &
    #          df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    # {
    #   vec.final.sequence <- c( vec.final.sequence,
    #                            "I" )
    # } # Add Janssen mutation
    # else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 658 &
    #          df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    # {
    #   vec.final.sequence <- c( vec.final.sequence,
    #                            "V" )
    # } 
    # # Add 368R mutation (Mangai's request)
    # else if( df.no_aln.no_hxb2.insert_hxb2.repaired_aa$no_hxb2[ i ] == 368 &
    #          df.no_aln.no_hxb2.insert_hxb2.repaired_aa$insert_hxb2[ i ] == " " )
    # {
    #   vec.final.sequence <- c( vec.final.sequence,
    #                            "R" )
    # }
    else
    {
      vec.final.sequence <- c( vec.final.sequence,
                               df.no_aln.no_hxb2.insert_hxb2.repaired_aa$aa_newSequence[ i ] )
    }
  }
  
  # Add N- and C-terminal tags
  vec.final.sequence.Nt.Ct <- c( unlist( strsplit( "MPMGSLQPLATLYLLGMLVASVLAQQAENL",
                                                   "" ) ),
                                 vec.final.sequence,
                                 unlist( strsplit( "GSAPTKAKRRVVQREKR",
                                                   "" ) ) )
  # vec.final.sequence.Nt.Ct <- c( vec.final.sequence )
  
  #==================================================
  # Save repaired sequence with additional mutations
  vec.final.sequence.Nt.Ct <- vec.final.sequence.Nt.Ct[ vec.final.sequence.Nt.Ct != "-" ]
  write.fasta( seqs = vec.final.sequence.Nt.Ct,
               id = output.prefix,
               file = paste( output.prefix,
                             "_DS-SOSIP.fasta",
                             sep = "" ) )
  
}



#==================================================
#==================================================
#==================================================
# Main

#==================================================
# start time measurement
start.main <- proc.time()
#==================================================

#==================================================
# Command line arguments
file.repair.sequence <- commandArgs()[ 3 ]
subtype <- commandArgs()[ 4 ]
exe.muscle <- commandArgs()[ 5 ]

#==================================================
# Source function

#==================================================
# Set working directory

#==================================================
# Set variables
file.reference.alignment <- "HIV1_FLT_2016_env_PRO_cleanLastCol_final.fasta"




#==================================================
# Run ADROITrimer function
ADROITrimer( file.reference.alignment,
             file.repair.sequence,
             subtype,
             exe.muscle )

  
  


#==================================================
# Measure and print out time needed for the script
end.main <- proc.time()
duration.main <- end.main-start.main
print( paste( "Script duration:", round( duration.main[3] / 60, 2 ), "min") )
#==================================================

# Main end
#==================================================
#==================================================
#==================================================