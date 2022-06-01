#!/usr/bin/perl
# #!/diskmnt/Software/LabCode/lyao/Tools/anaconda3/bin/perl
#
#  THIS SCRIPT IS FOR FINDING TISSUE-SPECIFIC EXPRESSION OF GENES FROM THE
#  GTEX DATABASE (USED FOR THE 2019 MM PAPER GENES IN FIG 7A)
#
#  ######################
#  #  PROGRAMMER NOTES  #
#  ######################
#
#  -------
#  GENERAL
#  -------
#
#     The analysis for this particular case was for a specific gene list
#     from Ruiyang. Therefore, I pre-processed the full GTex files to parse-out
#     only those genes of interest using a few 1-off throwaway scripts.
#     If you run this on the full GTex files, it will analyze *all* the genes
#     in there.
#
#  -----------------------
#  FILES AND THEIR FORMATS
#  -----------------------
#
#     EXPRESSION DATA FILE: columns and rows are GTex samples and gene names,
#     respectively, i.e. each column gives the expression of every gene for
#     that sample in units of Transcripts Per Million (TPM). All data are
#     tab-separated. For example:
#
#     Name  Description  GTEX-1117F-0226-SM-5GZZ7  GTEX-111CU-1826-SM-5GZYN ...
#     ENSG00000169504.10      CLIC4   174.7   162     270.2 ...
#     ENSG00000116209.7       TMEM59  141.7   139.7   112.4 ...
#     ENSG00000143575.10      HAX1    87.32   101.5   94.86 ...
#             :                :        :       :       :
#
#     META-DATA FILE: this is the file that links the GTex file name to its
#     tissue type. The parser here expects this file in the form of 2
#     tab-separated columns.
#
#     GTEX-1117F-0003-SM-58Q7G        Whole Blood
#     GTEX-1117F-0003-SM-5DWSB        Whole Blood
#     GTEX-1117F-0003-SM-6WBT7        Whole Blood
#     GTEX-1117F-0226-SM-5GZZ7        Adipose - Subcutaneous
#     GTEX-1117F-0426-SM-5EGHI        Muscle - Skeletal
#     GTEX-1117F-0526-SM-5EGHJ        Artery - Tibial
#               :                          :
#
#  -----------------
#  MAPPING TEST DATA
#  -----------------
#
#     These are a few test cases to make sure the mapping works correctly, i.e.
#     mapping column number in the expression file to the GTex sample name to
#     the tissue type
#
#     2 => GTEX-111FC-0226-SM-5N9B8 => Adipose
#     15 => GTEX-11EQ8-0226-SM-5EQ5G => Adipose
#     11684 => GTEX-ZVZP-0006-SM-51MSW => Blood
#     
############
#  SET UP  #
############

#__STANDARD PRAGMAS
   use v5.18.2;
   use strict;
   use warnings;

my ($outdir) = @ARGV;

if (not defined $outdir) {
  die "Need path of GTEX expression of DEGs from output dir \n";
}

#__STANDARD PERL/C-PAN LIBRARIES
   use Statistics::Distributions;
   use Math::SigFigs;

#__SPECIAL LIBRARIES
   use lib "./GTEX_lib/";
   #my $libdir;
   #BEGIN { $libdir = '../Scripts' }
   #use lib $libdir;
   use Statistics::Outliers;
###################
#  CONFIGURATION  #
###################

#__DIRECTORY WHERE THE PRE-PROCESSED GTEX SOURCE FILES LIVE
   #my $outdir = "/diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Results/V8/Union_discovery";

#__GTEX DATA FILE
   my $gtex_data_file = join ('/', $outdir, "GTEX_expression_matrix_subset_tumor_cells_DEG.txt");
#  my $gtex_data_file = "expression_data_GTEX_subset_Melanoma_DE_genes.txt";

#__GTEX ANNOTATION (META-DATA) FILE
   my $gtex_annot_file = join ('/', $outdir, "GTEX_metadata_subset_tumor_cells_DEG.txt");
#  my $gtex_annot_file = "metadata_GTEX_subset_Melanoma_DE_genes.txt";

#__SET GLOBAL FORMATTING VARIABLES FOR OUTPUT
   $^ = "OUTPUT_TOP";
   $~ = "OUTPUT";
   $- = 0; # closes out any previous report so current "TOP" will be used

####################
#  PRE-PROCESSING  #
####################

#__GET MAPPING: LIST POSITION TO SAMPLE
   my ($map_position_2_sample, $count_positions) =
         _mapping_position_2_sample_ ($gtex_data_file);

#__GET MAPPING: SAMPLE NAME TO TISSUE
   my ($map_sample_2_tissue, $count_samples) =
         _mapping_sample_2_tissue_ ($gtex_annot_file);

#__TEST CASES TO CONFIRM MAPPING IS WORKING CORRECTLY
#  for my $i (2, 15, 11684) {
#     my $sample = $map_position_2_sample->{$i};
#     my $tissue = $map_sample_2_tissue->{$sample};
#     print "$i   $sample   $tissue\n";
#  }

#__READ EXPRESSION DATA
   my $data = {};
   open (F, $gtex_data_file) || die "cant open file '$gtex_data_file'";
   while (<F>) {
      chomp;
      next if /^#/ || /^Name/;

   #__PARSE LINE AND SHIFT OFF 1ST ELEMENT AS JUNK AND 2ND AS GENE
      my @fields = split /\t/;
      shift @fields;
      my $gene = shift @fields;

   #__GO THROUGH THE DATA FIELD ULTIMATELY FINDING TISSUE FOR EACH EXPRESSION
      for (my $i = 0; $i < $#fields; $i++) {
         my $sample = $map_position_2_sample->{$i};
         my $tissue = $map_sample_2_tissue->{$sample};

      #__RECORD THIS DATUM
         push @{$data->{$gene}->{$tissue}}, $fields[$i];
      }
   }
   close (F);

###############
#  MAIN CODE  #
###############

#__DISCERN GOOD TEST CANDIDATES USING AN OUTLIER HEURISTIC
   my ($number_of_outliers, $total_possible_tests) = (0, 0);
   my ($hypothesis_tests, $background) = ({}, {});
   foreach my $gene (keys %{$data}) {

   #__LOOK AT EACH TISSUE TYPE IN THIS GENE
      my $averages_list = [];
      my $averages_hash = {};
      foreach my $tissue (keys %{$data->{$gene}}) {
         $total_possible_tests++;

      #__THE AVERAGE EXPRESSION VALUE FOR THIS TISSUE
         my $avg = _mean_ ($data->{$gene}->{$tissue});

      #__UNLIKELY THAT THERE ARE REPEATED AVERAGES SO USE THAT AS KEY
         if (defined $averages_hash->{$avg}) {
#            die "$gene: repeated avg for $tissue and $averages_hash->{$avg}";
            print "$gene: repeated avg for $tissue and $averages_hash->{$avg}"; 
         } else {
            $averages_hash->{$avg} = $tissue;
         }

      #__ADD TO LIST OF AVERAGES
         push @{$averages_list}, $avg;
      }

   #__GET HIGH-END OUTLIERS AS BEING THE BEST TEST CANDIDATES USING A HEURISTIC
      my ($low_rejects, $hi_rejects);
      ($averages_list, $low_rejects, $hi_rejects) =

         #  A LOT OF METHODS INCLUDING INTERQUARTILE RETURNED JUNK
         #
         #__NONE OF THESE WORKED WELL: MOSTLY LET THROUGH TOO MUCH JUNK
         #  Statistics::Outliers::tukey_iqr ($averages_list, 1); # PERMISSIVE
         #  Statistics::Outliers::tukey_iqr ($averages_list); # STANDARD
         #  Statistics::Outliers::median_absolute_dev ($averages_list); # BAD
         #  Statistics::Outliers::tukey_iqr ($averages_list, 2); # STRICTER
         #
         #__STRICTER BUT STILL ALLOWS SPLEEN, SALIVARY & INTESTINE FOR TNFRSF17
         #  Statistics::Outliers::tukey_iqr ($averages_list, 3);

         #__THIS GETS IT RIGHT FOR TNFRSF17: *** USE THIS ***
            Statistics::Outliers::peirce_criterion ($averages_list);

   #__MATCH THESE NON-OUTLIER AVERAGES BACK TO THEIR TISSUE TYPES FOR BACKGROUND
      foreach my $avg (@{$averages_list}) {
         my $tissue = $averages_hash->{$avg};
         push @{$background->{$gene}}, $tissue;
      }
      foreach my $avg (@{$low_rejects}) {
         my $tissue = $averages_hash->{$avg};
         push @{$background->{$gene}}, $tissue;
      }

   #__MATCH THESE HIGH-END AVERAGES BACK TO THEIR TISSUE TYPES FOR TESTING
      foreach my $avg (@{$hi_rejects}) {
         my $tissue = $averages_hash->{$avg};
         push @{$hypothesis_tests->{$gene}}, $tissue;
      }

   #------ FOR DEBUGGING
   # if ($gene eq "CD38") {
   #    print "BASIC BACKGROUND: ", join (' ', @{$background->{$gene}}), "\n";
   #    print "TEST TISSUE SET: ", join (' ', @{$hypothesis_tests->{$gene}}), "\n";
   # }
   #------ FOR DEBUGGING
   }

#__DIFFERENCE-OF-MEANS TEST FOR EACH MARKED COMBINATION OF TISSUE TYPE AND GENE
   my ($num_tests, $pval_results) = (0, {});
   foreach my $gene (keys %{$hypothesis_tests}) {

   #__TEST EACH OUTLIER CANDIDATE
      foreach my $tissue (@{$hypothesis_tests->{$gene}}) {

      #__BACKGROUND CONTAINS ALL EXPRESSIONS OF ALL NON-OUTLIER TISSUE TYPES
         # my ($num_bgrnd, $avg_bgrnd, $variance_bgrnd) = (0, 0, 0);
         my $express_bgrnd = [];
         foreach my $bgrnd_tissue (@{$background->{$gene}}) {
            push @{$express_bgrnd}, @{$data->{$gene}->{$bgrnd_tissue}};
         }

      #__BACKGROUND ALSO CONTAINS ALL OF THE OTHER OUTLIER TISSUE TYPES
         foreach my $other_tissue (@{$hypothesis_tests->{$gene}}) {
            next if $other_tissue eq $tissue;
            push @{$express_bgrnd}, @{$data->{$gene}->{$other_tissue}};
         }

      #__DESCRIPTIVE STATS OF THE OVERALL BACKGROUND
         my $num_bgrnd = scalar @{$express_bgrnd};
         my ($avg_bgrnd, $variance_bgrnd) = _mean_and_variance_($express_bgrnd);

      #__DESCRIPTIVE STATS OF THE TEST TISSUE
         my $num_test = scalar @{$data->{$gene}->{$tissue}};
         my ($avg_test, $variance_test) = _mean_and_variance_ (
            $data->{$gene}->{$tissue}
         );

      #__T-STATISTIC: NULL = TEST MEAN IS *NOT* HIGHER THAN BACKGROUND MEAN
         my ($t, $fold_change);
         if ($avg_test >= $avg_bgrnd) {
            $fold_change = $avg_test / $avg_bgrnd;
            $t = $avg_test - $avg_bgrnd;
            $t /= sqrt ($variance_test/$num_test + $variance_bgrnd/$num_bgrnd);
         } else {
            print "    $gene   $tissue   TEST ($avg_test) < BACKGROUND " .
                  "($avg_bgrnd) : NO TEST PERFORMED\n";
            next;
         }

      #__DEGREES OF FREEDOM
         my $dof = $num_bgrnd + $num_test - 2;

      #__HYPOTHESIS TEST
         my $p_value = Statistics::Distributions::tprob ($dof, $t);
#?         print "$gene  $tissue        T = $t   DOF = $dof   P = $p_value\n";

      #__SAVE TEST RESULTS
         push @{$pval_results->{$p_value}}, [$gene, $tissue, $fold_change];
         $num_tests++;
      }
   }

#__MULTIPLE TEST CORRECTION
   my ($output_by_fdr, $output_by_gene) = ([], {});
   my ($rank, $fdr_previous, $min_non_zero, $fdr_non_zero) = (1, 0, 99999, 0);
   foreach my $p_value (sort _numerical_ keys %{$pval_results}) {
      foreach my $list_reference (@{$pval_results->{$p_value}}) {
         my ($gene, $tissue, $fold_change) = @{$list_reference};

      #__ACTUAL FDR CALCULATION
         my $fdr = $p_value * $num_tests / $rank;

      #__MODIFIER 1: FDR DOES NOT EXCEED UNITY
         $fdr = 1 if $fdr > 1;

      #__MODIFIER 2: FDR IS MONOTONIC
         $fdr = $fdr_previous if $fdr < $fdr_previous;
         $fdr_previous = $fdr;

      #__RIGMAROLE ROUNDING FOR 2 DIGITS BEHIND THE DECIMAL IN SCI-NOTATION
         my $fdr_round = FormatSigFigs ($fdr, 5);
         my $fdr_sci_notation = sprintf ("%.4e", $fdr_round);

      #__RIGMAROLE ROUNDING FOR 2 DIGITS BEHIND THE DECIMAL IN SCI-NOTATION
         my $p_value_round = FormatSigFigs ($p_value, 5);
         my $p_value_sci_notation = sprintf ("%.4e", $p_value_round);

      #__PULL OUT THE SMALLEST NON-ZERO VALUE
         if ($fdr) {
            if ($fdr < $min_non_zero) {
               $min_non_zero = $fdr;
               $fdr_non_zero = $fdr_sci_notation;
            }
         }

      #__SAVE KEYED BY GENE
         push @{$output_by_gene->{$gene}},
            [$tissue, $fold_change, $p_value_sci_notation, $fdr_sci_notation];

      #__SAVE KEYED BY FDR
         if ($fdr) {
            push @{$list_reference}, $fdr_sci_notation;
         } else {
            push @{$list_reference}, $fdr;
         }
         push @{$output_by_fdr}, $list_reference;

      #__OUTPUT AND INCREMENT
#?         print "RANK $rank:  $list_reference->[0]   $list_reference->[1]    " .
#?               "PVAL = $p_value   FDR = $fdr_sci_notation\n";
#?         $rank++;
      }
   }

#__FINAL OUTPUT
   foreach my $gene (sort keys %{$output_by_gene}) {
      foreach my $tissue_datum (@{$output_by_gene->{$gene}}) {
         my ($tissue, $fold_change, $pval, $fdr) = @{$tissue_datum};
         my $fdr_print;
         if ($fdr) {
            $fdr_print = $fdr;
         } else {
            $fdr_print = "< $fdr_non_zero";
         }
#__FORMATTING
format OUTPUT_TOP =
  ----------------------------------------------------------------------------
  GENE           TISSUE          P-VALUE        FDR            FOLD DIFFERENCE
  ----------------------------------------------------------------------------
.
format OUTPUT =
  @<<<<<<<<<<<   @<<<<<<<<<<<<<  @<<<<<<<<<<<<  @<<<<<<<<<<<<  @<<<<<<
  $gene, $tissue, $pval, $fdr_print, $fold_change
.
         write;
      }
   }

#####################
#  POST-PROCESSING  #
#####################

#  NONE

################################################################################
#                                                                              #
#                                  SUBROUTINES                                 #
#                                                                              #
################################################################################

#__STANDARD PERL NUMERICAL SORT
   sub _numerical_ {$a <=> $b}

#############################
#  SUROUTINES:  STATISTICS  #
#############################

#  ----------------------
#  SIMPLE ARITHMETIC MEAN
#  ----------------------

sub _mean_ {
   my ($list) = @_;
   my $mean = 0;
   my $num_vals = scalar @{$list};
   foreach my $number (@{$list}) {
      $mean += $number;
   }
   $mean /= $num_vals;
   return $mean;
}

#  ---------------------------
#  MEAN AND VARIANCE
#  ---------------------------

sub _mean_and_variance_ {
   my ($list) = @_;

#__CENTRAL LOCATION IS THE MEAN
   my $x_1_bar = _mean_ ($list);

#__SUM OF SQUARES OF DIFFEFERENCES FROM MEAN
   my $sum_sq_of_diffs_1 = _sum_square_of_diffs_ ($list, $x_1_bar);

#__UNBIASED ESTIMATOR IS SIZE MINUS 1
   my $denominator = scalar @{$list} - 1;

#__SAMPLE VARIANCE
   my $variance;
   if ($denominator) {
      $variance = $sum_sq_of_diffs_1 / $denominator;
   } else {
      $variance = 0;
   }

#__RETURN MEAN AND STANDARD DEVIATION
   return ($x_1_bar, $variance);
}

#  -----------------------------
#  SUM OF SQUARES OF DIFFERENCES
#  -----------------------------

sub _sum_square_of_diffs_ {
   my ($list, $central_location) = @_;
   my $sum_square_of_diffs = 0;
   foreach my $number (@{$list}) {
      my $sq_diff = ($number - $central_location)*($number - $central_location);
      $sum_square_of_diffs += $sq_diff;
   }
   return $sum_square_of_diffs;
}

##########################
#  SUROUTINES:  MAPPING  #
##########################

#  --------------------------------
#  MAPPING: LIST POSITION TO SAMPLE
#  --------------------------------

   sub _mapping_position_2_sample_ {
      my ($file) = @_;
      my ($count, $mapping) = (0, {});
      open (F, $file) || die "cant open file '$file'";
      while (<F>) {
         chomp;

      #__PROCESS THE HEADER FROM WHICH WE CAN INFER THE MAPPING
         if (/^Name/) {

         #__PARSE LINE AND SHIFT OFF FIRST 2 ELEMENTS WHICH ARE NOT SAMPLES
            my @fields = split /\t/;
            shift @fields;
            shift @fields;

         #__DIAGNOSTICS
            my $count = scalar @fields;
#?            print "mapping_position_2_sample: processed $count samples\n";

         #__MAP LIST POSITION TO SAMPLE NAME
            for (my $i = 0; $i < $#fields; $i++) {
               $mapping->{$i} = $fields[$i];
            }

         #__DONE WITH FILE
            last;
         }
      }
      close (F);
      return ($mapping, $count);
   }

#  --------------------------------
#  MAPPING: SAMPLE TO TISSUE TYPE
#  --------------------------------
 
   sub _mapping_sample_2_tissue_ {
      my ($file) = @_;
      my ($count, $mapping) = (0, {});

   #__GET MAPPING FROM SPECIFIC TO GENERAL TISSUE TYPE
      my $conversion = &mapping_specific_2_general;

   #__READ META-DATA FILE
      open (F, $file) || die "cant open file '$file'";
      while (<F>) {

      #__READ A DATA LINE
         chomp;
         next if /^#/;
         my ($sample, $specific_tissue_type) = split /\t/;

      #__CONVERT THE SPECIFIC TISSUE TYPE TO HIGH-LEVEL GENERAL TYPE
         my $general_tissue_type;
         if (defined $conversion->{$specific_tissue_type}) {
            $general_tissue_type = $conversion->{$specific_tissue_type};
         } else {
            die "$specific_tissue_type is not a recognized tissue type";
         }

      #__MAPPING GENERAL TISSUE TYPE TO SAMPLE
         $mapping->{$sample} = $general_tissue_type;
         $count++;
      }
      close (F);
#?      print "mapping_sample_2_tissue: processed $count samples\n";
      return ($mapping, $count);
   }

#  -----------------------------------------------------
#  MAPPING: SPECIFIC TISSUE TYPE TO GENERAL TISSUE TYPE 
#  -----------------------------------------------------

   sub mapping_specific_2_general {
      my $exact_2_general_type = {
         'Adipose - Subcutaneous' => 'Adipose',
         'Adipose - Visceral (Omentum)' => 'Adipose',
         'Adrenal Gland' => 'Adrenal',
         'Artery - Aorta' => 'Artery',
         'Artery - Coronary' => 'Artery',
         'Artery - Tibial' => 'Artery',
         'Bladder' => 'Bladder',
         'Brain - Amygdala' => 'Brain',
         'Brain - Anterior cingulate cortex (BA24)' => 'Brain',
         'Brain - Caudate (basal ganglia)' => 'Brain',
         'Brain - Cerebellar Hemisphere' => 'Brain',
         'Brain - Cerebellum' => 'Brain',
         'Brain - Cortex' => 'Brain',
         'Brain - Frontal Cortex (BA9)' => 'Brain',
         'Brain - Hippocampus' => 'Brain',
         'Brain - Hypothalamus' => 'Brain',
         'Brain - Nucleus accumbens (basal ganglia)' => 'Brain',
         'Brain - Putamen (basal ganglia)' => 'Brain',
         'Brain - Spinal cord (cervical c-1)' => 'Brain',
         'Brain - Substantia nigra' => 'Brain',
         'Breast - Mammary Tissue' => 'Breast',
         'Cells - EBV-transformed lymphocytes' => 'lymphocytes',
         'Cells - Transformed fibroblasts' => 'Fibroblast',
         'Cervix - Ectocervix' => 'Vagina',
         'Cervix - Endocervix' => 'Vagina',
         'Colon - Sigmoid' => 'Colon',
         'Colon - Transverse' => 'Colon',
         'Esophagus - Gastroesophageal Junction' => 'Esophagus',
         'Esophagus - Mucosa' => 'Esophagus',
         'Esophagus - Muscularis' => 'Esophagus',
         'Fallopian Tube' => 'Ovary',
         'Heart - Atrial Appendage' => 'Heart',
         'Heart - Left Ventricle' => 'Heart',
         'Kidney - Cortex' => 'Kidney',
         'Liver' => 'Liver',
         'Lung' => 'Lung',
         'Minor Salivary Gland' => 'Salivary',
         'Muscle - Skeletal' => 'Muscle',
         'Nerve - Tibial' => 'Nerve',
         'Ovary' => 'Ovary',
         'Pancreas' => 'Pancreas',
         'Pituitary' => 'Pituitary',
         'Prostate' => 'Prostate',
         'Skin - Not Sun Exposed (Suprapubic)' => 'Skin',
         'Skin - Sun Exposed (Lower leg)' => 'Skin',
         'Small Intestine - Terminal Ileum' => 'Intestine',
         'Spleen' => 'Spleen',
         'Stomach' => 'Stomach',
         'Testis' => 'Testis',
         'Thyroid' => 'Thyroid',
         'Uterus' => 'Uterus',
         'Vagina' => 'Vagina',
         'Whole Blood' => 'Blood',
      };
      return $exact_2_general_type;
   }

################################################################################
#                                                                              #
#                                     DATA                                     #
#                                                                              #
################################################################################
