package Bio::IsotopeCorrection::Chemicals;
################################################################################
################################################################################
# Author:    Christian Jungreuthmayer
# Company:   Austrian Centre of Industrial Biotechnology (ACIB)
# Copyright: ACIB 2014, 2015
# License:   GNU General Public License (GPL)
#            Artistic License
################################################################################

use strict;
use warnings;

use Carp;
use Bio::IsotopeCorrection::NaturalIsotopes;

use constant CONSIDER_ZERO => 1e-8;

our $VERSION = '0.01';

my $g_cache_n_over_k;
my $g_cache_prob;

################################################################################
################################################################################
sub new
{
   my $whoami = _whoami();
   my $class  = shift;
   my $params = shift;
   my $self = {};
   $g_cache_n_over_k = {};
   $g_cache_prob     = {};


   die "ERROR: $whoami: It seems that new() was already called for this object" if ref $class;
   die "ERROR: $whoami: Provided parameter is not a hash reference" if defined $params && ref($params) ne "HASH";

   $self->{nat_ab_tracer} = $params->{nat_ab_tracer};

   ############################################################################
   ############################################################################
   if( defined $params->{nat_isotopes} )
   {
      if( ref $params->{nat_isotopes} ne 'Bio::IsotopeCorrection::NaturalIsotopes' )
      {
         die "ERROR: $whoami: provided object is not of type 'NaturalIsotopes'";
      }
      $self->{nat_isotopes} = $params->{nat_isotopes};
   }
   else
   {
      die "ERROR: $whoami: provide object of type 'NaturalIsotopes' with hash-parameter 'nat_isotopes'";
   }
   ############################################################################

   ############################################################################
   ############################################################################
   if( defined $params->{filename_chem} )
   {
      my $filename = $params->{filename_chem};
      warn "INFO: $whoami: going to open file '$filename' for reading ...";
      my $num_sets_external = 0;
      open my $fh, $filename or die "ERROR: couldn't open file '$filename' for reading: $!\n";

      _process_data_sets($fh, $filename, $self);

      close $fh;
   }
   else
   {
      die "ERROR: $whoami: provide name of data file containing information about precursor and fragement (hash-key: filename_chem)";
   }
   ############################################################################

   ############################################################################
   ############################################################################
   if( defined $params->{filename_tof} )
   {
      my $filename = $params->{filename_tof};
      warn "INFO: $whoami: going to open file '$filename' for reading ...";
      open my $fh, $filename or die "ERROR: couldn't open file '$filename' for reading: $!\n";
      _read_tof_values($fh, $filename, $self);
      close $fh;
   }
   else
   {
      die "ERROR: $whoami: provide name of data file containing measured TOF values to be corrected (hash-key: filename_tof)";
   }
   ############################################################################

   ############################################################################
   ############################################################################
   if( $params->{tracer_purity} )
   {
      $self->{tracer_purity} = 1;
      my $filename = $params->{tracer_purity};
      warn "INFO: $whoami: going to open file '$filename' for reading ...";
      open my $fh, $filename or die "ERROR: couldn't open file '$filename' for reading: $!\n";
      _read_purity_info($fh, $filename, $self);
      close $fh;
   }
   ############################################################################

   ############################################################################
   # do the hard work of the correction
   ############################################################################
   ($self->{c_iso_matrix}, $self->{tracer_Nn_pairs}) = _compute_tracer_iso_matrix( $self );
   _check_tof_values( $self );
   _compute_isotope_permutations( $self );
   $self->{sum_rel_masses} = _sum_up_masses( $self );
   $self->{compound_probabilities}  = _compute_compound_probabilities( $self );
   ############################################################################

   ############################################################################
   # expand isotope combinations if (im)purity of tracer needs to be
   # accounted for
   ############################################################################
   if( $self->{tracer_purity} )
   {
      _expand_combinations_by_purity( $self );
   }
   ############################################################################

   ############################################################################
   ############################################################################
   $self->{merged_hash} = _create_merged_hash( $self );
   ############################################################################

   $self = bless $self, $class;

   return $self;
}
################################################################################


################################################################################
################################################################################
sub _triangulate_matrix
{
   my $whoami   = _whoami();
   my $self     = shift;
   my $corr_mat = shift;
   my $tri_mat  = [];
   my $tri_meas = [];

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);

   my $tof_measured      = $self->{tof_measured};

   my $num_rows = @$corr_mat;
   my $num_cols = @{$corr_mat->[0]};
   die "ERROR: $whoami: matrix to be triangulated is not a square matrix\n" unless $num_rows == $num_cols;

   # create copy of correction matrix
   for( my $r = 0; $r < @$corr_mat; $r++ )
   {
      for( my $c = 0; $c < @{$corr_mat->[$r]}; $c++ )
      {
         $tri_mat->[$r][$c] = $corr_mat->[$r][$c];
      }
   }
   # create a copy of measured values
   for( my $r = 0; $r < @$tof_measured; $r++ )
   {
      $tri_meas->[$r] = $tof_measured->[$r]{val};
   }

   for( my $k = $num_rows - 1; $k > 0; $k-- )
   {
      die "ERROR: $whoami: element $k/$k of correction matrix was 0\n" if $tri_mat->[$k][$k] == 0;

      my $tmp1 = $tri_mat->[$k][$k];
      for( my $r = $k - 1; $r >= 0; $r-- )
      {
         my $tmp2 = $tri_mat->[$r][$k];
         for( my $c = $k; $c >= 0; $c-- )
         {
            $tri_mat->[$r][$c] -= $tri_mat->[$k][$c]*$tmp2/$tmp1;
         }
         $tri_meas->[$r] -= $tri_meas->[$k]*$tmp2/$tmp1;
      }
   }

   # _print_matrix( $tri_mat, "Triangular matrix:\n");
   # print "tri_measured: @$tri_meas\n";

   return $tri_mat, $tri_meas;
}
################################################################################


################################################################################
################################################################################
sub do_correction
{
   my $whoami = _whoami();
   my $self = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);

   my $tof_measured      = $self->{tof_measured};
   my $merged_hash       = $self->{merged_hash};

   ############################################################################
   # correction
   ############################################################################
   my $correction_value   = [];
   my $tof_corrected      = [];
   my $tof_corrected_nrm  = [];
   my $tof_corrected_hash = {};

   $self->{tof_corrected}      = $tof_corrected;
   $self->{tof_corrected_norm} = $tof_corrected_nrm;
   $self->{tof_corrected_hash} = $tof_corrected_hash;
   $self->{correction_value}   = $correction_value;

   # build correction matrix
   my $correction_matrix = [];
   my $index_map = {};
   for( my $m = 0; $m < @$tof_measured; $m++ )
   {
      my $N = $tof_measured->[$m]{N};
      my $n = $tof_measured->[$m]{n};
      my $key = $N . '_' . $n;

      $index_map->{$key} = $m;
      
      for( my $n = 0; $n <  @$tof_measured; $n++ )
      {
         $correction_matrix->[$m][$n] = 0.0;
      }
   }

   for( my $m = 0; $m < @$tof_measured; $m++ )
   {
      my $N = $tof_measured->[$m]{N};
      my $n = $tof_measured->[$m]{n};

      my $key = $N . '_' . $n;
      my $id_y = $index_map->{$key};

      # print "m=$m: N=$N, n=$n, key=$key id_$id_y\n";
      next unless exists $merged_hash->{$key};

      for( my $i = 0; $i < @{$merged_hash->{$key}}; $i++ )
      {
         my $index = $merged_hash->{$key}[$i];
         my ($tracer_N, $tracer_n) = _get_tracer_Nn_by_merged_index($self, $index);
         my $tracer_Nmass = $tracer_N*$self->{rel_mass_tracer_isotope};
         my $tracer_nmass = $tracer_n*$self->{rel_mass_tracer_isotope};
         my $key1 = $tracer_Nmass . '_' . $tracer_nmass;
         my $id_x = $index_map->{$key1};
         my $prob = $self->{compound_probabilities}[$index];

         $correction_matrix->[$id_y][$id_x] += $prob;
      }
   }

   _print_matrix( $correction_matrix, "Correction matrix:\n" );

   if( $self->{tracer_purity} )
   {
      my ($triCorrMatrix, $tof_tri) = $self->_triangulate_matrix( $correction_matrix);

      for( my $m = 0; $m < @$tof_tri; $m++ )
      {
         my $left = $tof_tri->[$m];
         for( my $n = 0; $n < $m; $n++ )
         {
            $left -= $triCorrMatrix->[$m][$n]*$tof_corrected->[$n];
         }
         $tof_corrected->[$m] = $left/$triCorrMatrix->[$m][$m];
         $tof_corrected->[$m] = 0 if $tof_corrected->[$m] < 0;
      }
      _print_matrix( $triCorrMatrix, "Correction matrix (triangular):\n" );
   }
   else
   {
      for( my $m = 0; $m < @$tof_measured; $m++ )
      {
         my $left = $tof_measured->[$m]{val};
         for( my $n = 0; $n < $m; $n++ )
         {
            $left -= $correction_matrix->[$m][$n]*$tof_corrected->[$n];
         }
         $tof_corrected->[$m] = $left/$correction_matrix->[$m][$m];
         $tof_corrected->[$m] = 0 if $tof_corrected->[$m] < 0;
      }
   }

   ############################################################################
   # normalize corrected areas to Icor[0] = Imeas[0]
   ############################################################################
   # if( $tof_measured->[0]{val} > 0 )
   if( $tof_corrected->[0] > 0 )
   {
      my $tmp_tof_corrected = $tof_corrected->[0];
      for( my $m = 0; $m < @$tof_corrected; $m++ )
      {
         # print "   m=$m tof_corrected->[$m]=$tof_corrected->[$m]\n";
         $tof_corrected->[$m] = $tof_corrected->[$m]/ $tmp_tof_corrected*$tof_measured->[0]{val};
         # print "   m=$m tof_corrected->[$m]=$tof_corrected->[$m]\n";
      }
   }
   ############################################################################

   ############################################################################
   # normalize corrected areas
   ############################################################################
   my $corr_area_sum = 0;
   for( my $m = 0; $m < @$tof_corrected; $m++ )
   {
      $corr_area_sum += $tof_corrected->[$m];
   }
   $corr_area_sum = 1 if $corr_area_sum == 0;
   for( my $m = 0; $m < @$tof_corrected; $m++ )
   {
      $tof_corrected_nrm->[$m] = $tof_corrected->[$m]/$corr_area_sum;
   }
   ############################################################################

   ############################################################################
   # compute enrichment
   ############################################################################
   my $enrichment = 0;
   if( $self->{isotopologue} )
   {
      for( my $m = 1; $m < @$tof_corrected_nrm; $m++ )
      {
         $enrichment += $tof_corrected_nrm->[$m]*$m;
      }
      $enrichment /= scalar(@$tof_corrected_nrm) - 1;
   }
   else
   {
      # note: the computation of an enrichment for non-isotopologues
      # (precursor and fragment are different) does not really make
      # sens -> skipped
   }

   $self->{enrichment} = $enrichment;
   ############################################################################
}
################################################################################


################################################################################
################################################################################
sub _compute_correction_value
{
   my $whoami = _whoami();
   my $self  = shift;
   my $N     = shift;
   my $n     = shift;
   my $cor_val = 0;

   my $merged_hash        = $self->{merged_hash};
   my $tof_measured_hash  = $self->{tof_hash};
   my $tof_corrected_hash = $self->{tof_corrected_hash};

   my $key = $N . '_' . $n;

   ############################################################################
   ############################################################################
   unless( $tof_measured_hash->{$key} )
   {
      die "ERROR: $whoami: the combination (N=$N,n=$n) was not found in list of measured TOF values\n";
   }
   ############################################################################

   ############################################################################
   # check if we do need to correct this value in the first place
   ############################################################################
   return 0 unless exists $merged_hash->{$key};
   ############################################################################

   ############################################################################
   ############################################################################
   for( my $i = 0; $i < @{$merged_hash->{$key}}; $i++ )
   {
      my $index = $merged_hash->{$key}[$i];
      my ($tracer_N, $tracer_n) = _get_tracer_Nn_by_merged_index($self, $index);
      # print "DEBUG: $whoami: i=$i, k=$key, index=$index, tracer_N=$tracer_N, tracer_n=$tracer_n\n";

      my $tracer_Nmass = $tracer_N*$self->{rel_mass_tracer_isotope};
      my $tracer_nmass = $tracer_n*$self->{rel_mass_tracer_isotope};

      my $correct_tof = $tof_corrected_hash->{$tracer_Nmass . '_' . $tracer_nmass};

      # die "$whoami: i=$i: couldn't find value in hash 'tof_corrected_hash' for N=$tracer_N n=$tracer_n\n" unless defined $correct_tof;
      next unless defined $correct_tof;

      my $prob = $self->{compound_probabilities}[$index];

      print "DEBUG: $whoami: N=$N, n=$n: index=$index, tracer_N=$tracer_N(mass=$tracer_Nmass) tracer_n=$tracer_n(mass=$tracer_nmass) correct_tof=$correct_tof prob=$prob\n";

      $cor_val -= $prob*$correct_tof;
   }
   ############################################################################

   return $cor_val;
}
################################################################################


################################################################################
################################################################################
sub print_enrichment
{
   my $whoami = _whoami();
   my $self = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);

   if( $self->is_isotopologue() )
   {
      print "mean enrichment: " . $self->get_enrichment() . "\n";
   }
   else
   {
      print "WARNING: we are not analyzing an isotopologue -> concept of enrichment does not make sense!\n";
   }
}
################################################################################


################################################################################
################################################################################
sub get_precursor_as_string
{
   my $whoami     = _whoami();
   my $self       = shift;
   my $pre_string = '';

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);

   $pre_string .= $self->{tracer_elem} . $self->{precursor}{tracer};
 
   foreach my $elem ( sort keys(%{$self->{precursor}{elements}}) )
   {  
      $pre_string .= $elem . $self->{precursor}{elements}{$elem};
   }
 
   foreach my $elem ( sort keys(%{$self->{excluded_precursor}{elements}}) )
   {  
      $pre_string .= $elem . $self->{excluded_precursor}{elements}{$elem};
   }

   return $pre_string
}
################################################################################


################################################################################
################################################################################
sub get_fragment_as_string
{
   my $whoami     = _whoami();
   my $self       = shift;
   my $fra_string = '';

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);

   $fra_string .= $self->{tracer_elem} . $self->{fragment}{tracer};
 
   foreach my $elem ( sort keys(%{$self->{fragment}{elements}}) )
   {  
      $fra_string .= $elem . $self->{fragment}{elements}{$elem};
   }
 
   foreach my $elem ( sort keys(%{$self->{excluded_fragment}{elements}}) )
   {  
      $fra_string .= $elem . $self->{excluded_fragment}{elements}{$elem};
   }

   return $fra_string
}
################################################################################


################################################################################
################################################################################
sub print_result
{
   my $whoami = _whoami();
   my $self = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);

   my $tof_measured      = $self->{tof_measured};
   my $correction_value  = $self->{correction_value};
   my $tof_corrected     = $self->{tof_corrected};
   my $tof_corrected_nrm = $self->{tof_corrected_norm};

   my $pre_string = $self->get_precursor_as_string();
   my $fra_string = $self->get_fragment_as_string();

   print "Tracer-Isotope: ", $self->{tracer_iso}, "\n";
   print "Precursor: ", $pre_string, "\n";
   print "Fragment:  ", $fra_string, "\n";
   for( my $m = 0; $m <  @$tof_measured; $m++ )
   {
      my $N = $tof_measured->[$m]{N};
      my $n = $tof_measured->[$m]{n};
      # print "$m: (M=$N,m=$n): measured=$tof_measured->[$m]{val} corrected=$tof_corrected->[$m] correct_normalized=$tof_corrected_nrm->[$m] correction_value=$correction_value->[$m]\n";
      print "$m: (M=$N,m=$n): measured=$tof_measured->[$m]{val} corrected=$tof_corrected->[$m] correct_normalized=$tof_corrected_nrm->[$m]\n";
   }
}
################################################################################


################################################################################
################################################################################
sub is_isotopologue
{
   my $whoami = _whoami();
   my $self = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);

   return $self->{isotopologue};
}
################################################################################


################################################################################
################################################################################
sub get_enrichment
{
   my $whoami = _whoami();
   my $self = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);

   return $self->{enrichment};
}
################################################################################


################################################################################
################################################################################
sub get_correct_areas
{
   my $whoami = _whoami();
   my $self = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);

   return $self->{tof_corrected};
}
################################################################################


################################################################################
################################################################################
sub get_correct_areas_normalized
{
   my $whoami = _whoami();
   my $self = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);

   return $self->{tof_corrected};
}
################################################################################


################################################################################
################################################################################
sub _compute_isotope_permutations
{
   my $whoami       = _whoami();
   my $self         = shift;
   my $nat_isotopes = $self->{nat_isotopes};

   print "INFO: $whoami: entered\n";

   ############################################################################
   # loop over all elements
   ############################################################################
   my $arr_merged_isos;
   my $arr_merged_wghts;
   my $arr_iso_names;
   my $arr_iso_rel_masses;
   my $arr_elem_names;
   foreach my $elem (sort keys %{$self->{precursor}{elements}})
   {
      my $arr_Nn_pairs;
      print "INFO: $whoami: element '$elem'\n";
      my $sum_N = $self->{precursor}{elements}{$elem};
      my $sum_n = $self->{fragment}{elements}{$elem};
      print "INFO: $whoami: sum_N=$sum_N, sum_n=$sum_n\n";

      my $lowest_isotope = $nat_isotopes->get_lowest_isotope($elem);
      print "INFO: $whoami: lowest_isotope isotope '$lowest_isotope'\n";

      foreach my $iso (sort $nat_isotopes->get_all_isotopes_for_element($elem))
      {
         next if $iso eq $lowest_isotope;
         my $rel_mass_number = $nat_isotopes->get_massnumber_for_isotope_relative_to_lowest_isotope($iso);

         my $Nn_pairs = _compute_Nn_pairs($self, $elem, $rel_mass_number);

         push @$arr_Nn_pairs, $Nn_pairs;
         push @$arr_iso_names, $iso;
         push @$arr_iso_rel_masses, $rel_mass_number;
      }

      my ($merged_isos, $merged_wghts) = _merge_iso_pairs( $self, $elem, $arr_Nn_pairs, $sum_N, $sum_n);
      push @$arr_merged_isos, $merged_isos;
      push @$arr_merged_wghts, $merged_wghts;
      push @$arr_elem_names, $elem;
   }
   ############################################################################

   ############################################################################
   ############################################################################
   if( $self->{nat_ab_tracer} )
   {
      my $arr_Nn_pairs;
      # my $elem = 'C'; # we currently only support C as tracer material
      my $elem = $self->{tracer_elem};
      print "INFO: $whoami: element '$elem' for considering natural abundance of tracer\n";
      my $sum_N = $self->{precursor}{tracer};
      my $sum_n = $self->{fragment}{tracer};
      print "INFO: $whoami: sum_N=$sum_N, sum_n=$sum_n\n";

      my $lowest_isotope = $nat_isotopes->get_lowest_isotope($elem);
      print "INFO: $whoami: lowest_isotope isotope '$lowest_isotope'\n";

      foreach my $iso (sort $nat_isotopes->get_all_isotopes_for_element($elem))
      {
         next if $iso eq $lowest_isotope;
         my $rel_mass_number = $nat_isotopes->get_massnumber_for_isotope_relative_to_lowest_isotope($iso);

         my $Nn_pairs = _compute_Nn_pairs_nat_ab_tracer($self, $elem, $rel_mass_number);

         push @$arr_Nn_pairs, $Nn_pairs;
         push @$arr_iso_names, 'nat_ab_tracer_' . $iso;
         push @$arr_iso_rel_masses, $rel_mass_number;
      }

      my ($merged_isos, $merged_wghts) = _merge_iso_pairs( $self, $elem, $arr_Nn_pairs, $sum_N, $sum_n);
      push @$arr_merged_isos, $merged_isos;
      push @$arr_merged_wghts, $merged_wghts;
      push @$arr_elem_names, 'nat_ab_tracer';
   }
   ############################################################################

   ############################################################################
   ############################################################################
   my $arr_tracer;
   push @$arr_tracer, $self->{tracer_Nn_pairs};
   my ($merged_tracer, $merged_tracer_wghts) =  _merge_iso_pairs( $self, 'tracer', $arr_tracer, $self->{precursor}{tracer}, $self->{fragment}{tracer});
   push @$arr_merged_isos, $merged_tracer;
   push @$arr_merged_wghts, $merged_tracer_wghts;
   push @$arr_iso_names, 'tracer';
   push @$arr_elem_names, 'tracer';
   push @$arr_iso_rel_masses, $self->{rel_mass_tracer_isotope};
   ############################################################################


   ############################################################################
   ############################################################################
   my $merged_elems = _merge_elems($self, $arr_merged_isos, $arr_merged_wghts, $arr_elem_names);
   $self->{merged_elems}          = $merged_elems;
   $self->{merged_iso_names}      = $arr_iso_names;
   $self->{merged_iso_rel_masses} = $arr_iso_rel_masses;
   # print "$whoami: iso names: @$arr_iso_names\n";
   # print "$whoami: iso relatives masses: @$arr_iso_rel_masses\n";
   ############################################################################

   print "INFO: $whoami: leaving\n";
}
################################################################################


################################################################################
################################################################################
sub _create_merged_hash
{
   my $whoami = _whoami();
   my $self = shift;
   my $merged_hash;

   my $iso_names      = $self->{merged_iso_names};
   my $merged_elems   = $self->{merged_elems};
   my $sum_rel_masses = $self->{sum_rel_masses};
   my $tof_hash       = $self->{tof_hash};


   ############################################################################
   ############################################################################
   for( my $r = 0; $r < @$merged_elems; $r++ )
   {
      my $key = $sum_rel_masses->[$r][0] . '_' . $sum_rel_masses->[$r][1];

      #########################################################################
      # check if all natural isotopes are zero
      #########################################################################
      my $all_zero = 1;
      my $num_entries = @{$merged_elems->[$r]};
      for( my $e = 0; $e < $num_entries - 1; $e++ )
      {
         my $elem_name;
         if( $iso_names->[$e] =~ /^([A-Z][a-z]*)(\d+)$/ )
         {
            $elem_name = $1;
         }
         elsif( $iso_names->[$e] =~ /^nat_ab_tracer_([A-Z][a-z]*)(\d+)$/ )
         {
            $elem_name = $1;
         }
         elsif( $iso_names->[$e] =~ /^purity_([A-Z][a-z]*)(\d+)$/ )
         {
            $elem_name = $1;
         }
         else
         {
            die "ERROR: $whoami: r=$r, e=$e/$num_entries: couldn't retrieve element name for isotope '$iso_names->[$e]'\n";
         }
         
         my $N = $merged_elems->[$r][$e][0];
         my $n = $merged_elems->[$r][$e][1];

         # if( $N != 0 )
         # {
         #    $all_zero = 0;
         #    last;
         # }
      }
      #########################################################################

      #########################################################################
      #########################################################################
      # if( $all_zero == 0 )
      # {
         push @{$merged_hash->{$key}}, $r;
      # }
      #########################################################################
   }

   ############################################################################
   ############################################################################
   foreach my $key (sort keys %$merged_hash )
   {
      for( my $i = 0; $i < @{$merged_hash->{$key}}; $i++ )
      {
         my $index = $merged_hash->{$key}[$i];
         my ($tracer_N, $tracer_n) = _get_tracer_Nn_by_merged_index($self, $index);

         my $measured_tof = $tof_hash->{$tracer_N . '_' . $tracer_n};
      }
   }
   ############################################################################

   return $merged_hash;
}
################################################################################


################################################################################
################################################################################
sub _get_tracer_Nn_by_merged_index
{
   my $self  = shift;
   my $index = shift;

   my $merged_elems = $self->{merged_elems};
   my $num_entries  = @{$merged_elems->[$index]};

   my $N = $merged_elems->[$index][$num_entries - 1][0];
   my $n = $merged_elems->[$index][$num_entries - 1][1];

   return $N, $n;
}
################################################################################


################################################################################
# compute the probability for each combination of isotopes
################################################################################
sub _compute_compound_probabilities
{
   my $whoami = _whoami();
   my $self = shift;
   my $compound_probabilites = [];
   my $rel_intensities = [];
   my $tracer_idx = 1;

   my $nat_isotopes = $self->{nat_isotopes};
   my $iso_names = $self->{merged_iso_names};
   my $merged_elems = $self->{merged_elems};

   # print "INFO: $whoami: entered\n";

   ############################################################################
   ############################################################################
   for( my $i = 0; $i < @$iso_names; $i++ )
   {
      # print "INFO: $whoami: $i: iso_name: $iso_names->[$i]\n";
      if( $iso_names->[$i] eq 'tracer')
      {
         $tracer_idx = $i;
         next;
      }

      my $rel_intensity;

      if( $iso_names->[$i] =~ /^nat_ab_tracer_([A-Z][a-z]?\d+)/ )
      {
         my $iso_name_nat_ab_tracer = $1;
         # $rel_intensity = $nat_isotopes->get_relative_intensity($iso_name_nat_ab_tracer)/100.0;
         $rel_intensity = $nat_isotopes->get_relative_intensity($iso_name_nat_ab_tracer);
      }
      else
      {
         # get relative intensity and convert it from percent to number
         $rel_intensity = $nat_isotopes->get_relative_intensity($iso_names->[$i]);
      }

      # print "$whoami: $i: iso_name: $iso_names->[$i], rel_intensity=$rel_intensity\n";
      push @$rel_intensities, $rel_intensity;
   }
   ############################################################################


   ############################################################################
   ############################################################################
   for( my $r = 0; $r < @$merged_elems; $r++ )
   {
      if( $r%10000 == 0 )
      {
         my $memGB = _get_used_memory();
         print "INFO: $whoami: r=$r/",scalar(@$merged_elems)," memory=${memGB}GB\n";
      }

      # print "INFO: $whoami: r=$r: processing compound:\n";
      # _print_compound($merged_elems->[$r], $iso_names);

      my $num_entries = @{$merged_elems->[$r]};
      my %tot_sum_precursor;
      my %tot_sum_fragment;
      my $hash_key = '';

      #########################################################################
      # loop over all entries to get total number of atoms (over all isotopes)
      # for the precursor and the fragement
      #########################################################################
      for( my $e = 0; $e < $num_entries - 1; $e++ )
      {
         my $elem_name;
         my $elem_name_hash;
         if( $iso_names->[$e] =~ /^([A-Z][a-z]?)(\d+)$/ )
         {
            $elem_name      = $1;
            $elem_name_hash = $1;
            
         }
         elsif( $iso_names->[$e] =~ /^nat_ab_tracer_([A-Z][a-z]?)(\d+)$/ )
         {
            $elem_name      = $1;
            $elem_name_hash = 'nat_ab_tracer_' . $1;
         }
         else
         {
            die "ERROR: $whoami: couldn't retrieve element name for isotope '$iso_names->[$e]'";
         }
         $tot_sum_precursor{$elem_name_hash} += $merged_elems->[$r][$e][0];
         $tot_sum_fragment{$elem_name_hash}  += $merged_elems->[$r][$e][1];
         $hash_key .= $merged_elems->[$r][$e][0] . 'x' . $merged_elems->[$r][$e][1] . '_';
         # print "e=$e: elem_name=$elem_name\n";
         # print "e=$e: elem_name_hash=$elem_name_hash\n";
      }
      # check if we have already computed the probablity for an identical data set
      # if so, use this value
      #if( $g_cache_prob->{$hash_key} )
      #{
      #   $compound_probabilites->[$r] = $g_cache_prob->{$hash_key};
      #   next;
      #}
      #########################################################################

      #########################################################################
      # compute probability
      #########################################################################
      my $prob  = 1.0;
      my $last_elem_name = "";
      my $sum_num_precursor  = 0;
      my $sum_num_fragment   = 0;
      my $taken_tot_precursor = -1;
      my $were_in_nat_ab_tracer = 0;
      my $am_in_nat_ab_tracer;

      for( my $e = 0; $e < $num_entries - 1; $e++ )
      {
         # print "INFO: $whoami: r=$r, e=$e, num_entries=$num_entries iso_name=$iso_names->[$e]\n";
         ######################################################################
         # get element name from isotope
         ######################################################################
         my $elem_name;
         my $elem_name_hash;
         if( $iso_names->[$e] =~ /^([A-Z][a-z]?)(\d+)$/ )
         {
            $elem_name           = $1;
            $elem_name_hash      = $1;
            $am_in_nat_ab_tracer = 0;
         }
         elsif( $iso_names->[$e] =~ /^nat_ab_tracer_([A-Z][a-z]?)(\d+)$/ )
         {
            $elem_name           = $1;
            $elem_name_hash      = 'nat_ab_tracer_' . $1;
            $am_in_nat_ab_tracer = 1;
         }
         else
         {
            die "ERROR: $whoami: couldn't retrieve element name for isotope '$iso_names->[$e]'\n";
         }
         # print "INFO: $whoami: elem_name=$elem_name\n";
         ######################################################################
         
         ######################################################################
         # set variables used in following section
         ######################################################################
         my $tot_precursor;
         my $tot_fragment;
         my $num_precursor;
         my $num_fragment;
         

         if( $iso_names->[$e] =~ /^nat_ab_tracer_/ )
         {
            $tot_precursor = $self->{precursor}{tracer} - $merged_elems->[$r][$tracer_idx][0];
            $tot_fragment  = $self->{fragment}{tracer}  - $merged_elems->[$r][$tracer_idx][1];
            $num_precursor = $merged_elems->[$r][$e][0];
            $num_fragment  = $merged_elems->[$r][$e][1];
            # print "   tracer_idx=$tracer_idx\n";
            # print "   self->{precursor}{tracer}:$self->{precursor}{tracer}\n";
            # print "   self->{fragment}{tracer}: $self->{fragment}{tracer}\n";
            # print "   tracerN=$merged_elems->[$r][$tracer_idx][0] tracern=$merged_elems->[$r][$tracer_idx][1]\n";
         }
         else
         {
            $tot_precursor = $self->{precursor}{elements}{$elem_name};
            $tot_fragment  = $self->{fragment}{elements}{$elem_name};
            $num_precursor = $merged_elems->[$r][$e][0];
            $num_fragment  = $merged_elems->[$r][$e][1];
         }
         # print "   tot_precursor=$tot_precursor tot_fragment=$tot_fragment\n";
         my $num_lowest_isotope_precursor = $tot_precursor - $tot_sum_precursor{$elem_name_hash};
         my $num_lowest_isotope_fragment  = $tot_fragment  - $tot_sum_fragment{$elem_name_hash};
         # print "   num_lowest_isotope_fragment=$num_lowest_isotope_fragment\n";
         # print "   num_lowest_isotope_precursor=$num_lowest_isotope_precursor\n";

         ######################################################################
         # if we deal with a new element reset some of the variables
         # and determine part of probability that is caused by most
         # abundant isotope which is not stored explicitely in arrays
         ######################################################################
         # print "   last_elem_name=$last_elem_name, elem_name=$elem_name\n";
         # print "   were_in_nat_ab_tracer=$were_in_nat_ab_tracer, am_in_nat_ab_tracer=$am_in_nat_ab_tracer\n";
         if( $last_elem_name ne $elem_name || $were_in_nat_ab_tracer != $am_in_nat_ab_tracer )
         {
            # print "   RESET counters\n";
            $sum_num_precursor = 0;
            $sum_num_fragment  = 0;

            my $lowest_isotope = $nat_isotopes->get_lowest_isotope($elem_name);
            my $rel_intensity_lowest_isotope = $nat_isotopes->get_relative_intensity($lowest_isotope);

            my $num_lowest_isotope = $tot_precursor - $tot_sum_precursor{$elem_name_hash};

            $taken_tot_precursor = $tot_precursor;

            $prob *= _n_over_k($tot_fragment, $num_lowest_isotope_fragment)*($rel_intensity_lowest_isotope**$num_lowest_isotope);
            # print "   rel_intensity_lowest_isotope=$rel_intensity_lowest_isotope\n";
            # print "   tot_fragment=$tot_fragment, num_lowest_isotope_fragment=$num_lowest_isotope_fragment num_lowest_isotope=$num_lowest_isotope\n";
            # print "   n_over_k($tot_fragment, $num_lowest_isotope_fragment)=",_n_over_k($tot_fragment, $num_lowest_isotope_fragment);
            # print "   prob=$prob\n";

            for( my $s = 0; $s < $num_lowest_isotope_fragment; $s++ )
            {
               $prob  *= ($num_lowest_isotope_precursor - $s)/$taken_tot_precursor;
               # print "      s=$s, prob=$prob\n";
               $taken_tot_precursor--;
            }
            $sum_num_fragment += $num_lowest_isotope_fragment;
         }
         # print "   AAA: prob=$prob\n";
         ######################################################################

         die "INTERNAL ERROR: $whoami: invalid value for taken_tot_precursor ($taken_tot_precursor)\n" if $taken_tot_precursor < 0;

         ######################################################################
         # determine probability that is caused by this isotope
         ######################################################################
         # print "   num_precursor=$num_precursor taken_tot_precursor=$taken_tot_precursor num_fragment=$num_fragment\n";
         for( my $s = 0; $s < $num_fragment; $s++ )
         {
            # print "      s=$s num_precursor=$num_precursor taken_tot_precursor=$taken_tot_precursor\n";
            $prob  *= ($num_precursor - $s)/$taken_tot_precursor;
            $taken_tot_precursor--;
         }
         # print "   BBB: prob=$prob\n";
         ######################################################################

         $tot_precursor -= $sum_num_precursor;
         $tot_fragment  -= $sum_num_fragment;
         # print "   sum_num_precursor=$sum_num_precursor sum_num_fragment=$sum_num_fragment\n";
         # print "   tot_precursor=$tot_precursor tot_fragment=$tot_fragment\n";

         ######################################################################
         ######################################################################
         my $N_over_K      = _n_over_k($tot_precursor,$num_precursor);
         my $n_over_k      = _n_over_k($tot_fragment,$num_fragment);
         # print "   N_over_K=$N_over_K, n_over_k=$n_over_k\n";
         $prob  *= $N_over_K*($rel_intensities->[$e]**$num_precursor)*$n_over_k;  # old feature
         # print "   prob=$prob\n";
         ######################################################################

         ######################################################################
         # update variable values for next iteration of this loop
         ######################################################################
         $last_elem_name     = $elem_name;
         $sum_num_precursor += $num_precursor;
         $sum_num_fragment  += $num_fragment;
         ######################################################################
         $were_in_nat_ab_tracer = $am_in_nat_ab_tracer;
      }


      #if( $self->{nat_ab_tracer} )
      #{
      #   # we consider natural abundance of not label fraction
      #}

      $compound_probabilites->[$r] = $prob;
      $g_cache_prob->{$hash_key}   = $prob;
      # print "$whoami: r=$r: prob=$prob\n";
      #########################################################################
   }

   # print "INFO: $whoami: leaving\n";

   return $compound_probabilites;
}
################################################################################


################################################################################
################################################################################
sub _faculty
{
   my $whoami = _whoami();
   my $inp = shift;

   die "$whoami: invalid input '$inp'\n" if $inp < 0;

   $inp == 0 ? return 1 : return $inp*_faculty($inp - 1);
}
################################################################################


################################################################################
################################################################################
sub _expand_combinations_by_purity
{
   my $whoami = _whoami();
   my $self = shift;

   my $merged_elems            = $self->{merged_elems};
   my $compound_probs          = $self->{compound_probabilities};
   my $sum_rel_masses          = $self->{sum_rel_masses};
   my $iso_names               = $self->{merged_iso_names};
   my $nat_isotopes            = $self->{nat_isotopes};
   my $rel_mass_tracer_isotope = $self->{rel_mass_tracer_isotope};

   my $num_iso_combinations = @$merged_elems;
   my $num_new_iso_combinations = 0;

   print "INFO: $whoami: expanding isotope combinations because of impurity of tracer ...\n";

   my @all_tracer_elem_isotopes = $nat_isotopes->get_all_isotopes_for_element($self->{tracer_elem});
   my @purities;
   for( my $i = @all_tracer_elem_isotopes - 1; $i >= 0; $i-- )
   {
      $purities[$i] = $self->{purity}{$all_tracer_elem_isotopes[$i]};
      unshift @$iso_names, 'purity_' . $all_tracer_elem_isotopes[$i];
   }

   for( my $r = 0; $r < $num_iso_combinations; $r++ )
   {
      my $num_entries = @{$merged_elems->[$r]};

      my $e_tracer = @{$merged_elems->[$r]} - 1;
      my $tracer_N = $merged_elems->[$r][$e_tracer][0];
      my $tracer_n = $merged_elems->[$r][$e_tracer][1];
      # my $tracer_N = $self->{precursor}{tracer};
      # my $tracer_n = $self->{fragment}{tracer};

      my ( $merged_purity, $merged_purity_wghts ) = _compute_combined_purity_Nn_pairs( $self, $tracer_N, $tracer_n );

      #########################################################################
      # new elements with values not equal to N=0 and n=0
      #########################################################################
      for( my $i = 1; $i < @$merged_purity; $i++ )
      {
         my $new_idx = $num_iso_combinations + $num_new_iso_combinations;

         ######################################################################
         # copy all existing entries to new merged elements
         ######################################################################
         for( my $e = 0; $e < $num_entries; $e++ )
         {
            $merged_elems->[$new_idx][$e][0]   = $merged_elems->[$r][$e][0];
            $merged_elems->[$new_idx][$e][1]   = $merged_elems->[$r][$e][1];
         }
         ######################################################################

         ######################################################################
         # now, unshift new entries because of tracer purity
         ######################################################################
         my $sum_N = 0;
         my $sum_n = 0;
         my $sum_N_wgt = 0;
         my $sum_n_wgt = 0;
         my $prop_corrector = 1;
         my $free_N_spots = $tracer_N;
         my $free_n_spots = $tracer_n;
         my $N_left       = $tracer_N;
         for( my $x =  @{$merged_purity->[$i]} - 1; $x >= 0; $x-- )
         {
            my $N = $merged_purity->[$i][$x][0];
            my $n = $merged_purity->[$i][$x][1];
            my $arr_ref;
            $arr_ref->[0] = $N;
            $arr_ref->[1] = $n;
            unshift @{$merged_elems->[$new_idx]}, $arr_ref;
            $sum_N += $N;
            $sum_n += $n;
            $prop_corrector *= $purities[$x]**$N*_n_over_k($free_N_spots, $N)*_n_over_k($free_n_spots, $n); # todo
            for( my $s = 0; $s < $n; $s++ )
            {
               $prop_corrector *= ($N - $s)/$N_left;
               $N_left--;
            }
            $free_N_spots -= $N;
            $free_n_spots -= $n;
         }

         $compound_probs->[$new_idx] = $compound_probs->[$r]*$prop_corrector;
         $sum_rel_masses->[$new_idx][0] = $sum_rel_masses->[$r][0] + $merged_purity_wghts->[$i][0] - $tracer_N*$rel_mass_tracer_isotope;
         $sum_rel_masses->[$new_idx][1] = $sum_rel_masses->[$r][1] + $merged_purity_wghts->[$i][1] - $tracer_n*$rel_mass_tracer_isotope;
         $num_new_iso_combinations++;
      }
      #########################################################################

      #########################################################################
      # now, update the original combiniation
      #########################################################################
      my $prop_corrector = 1;
      my $free_N_spots = $tracer_N;
      my $free_n_spots = $tracer_n;
      my $N_left       = $tracer_N;
      for( my $x =  @{$merged_purity->[0]} - 1; $x >= 0; $x-- )
      {
         my $N = $merged_purity->[0][$x][0];
         my $n = $merged_purity->[0][$x][1];
         my $arr_ref;
         $arr_ref->[0] = $N;
         $arr_ref->[1] = $n;
         unshift @{$merged_elems->[$r]}, $arr_ref;
         $prop_corrector *= $purities[$x]**$N*_n_over_k($free_N_spots, $N)*_n_over_k($free_n_spots, $n); # todo
         for( my $s = 0; $s < $n; $s++ )
         {
            $prop_corrector *= ($N - $s)/$N_left;
            $N_left--;
         }
         $free_N_spots -= $N;
         $free_n_spots -= $n;
      }
      $compound_probs->[$r] = $compound_probs->[$r]*$prop_corrector;
      $sum_rel_masses->[$r][0] = $sum_rel_masses->[$r][0] + $merged_purity_wghts->[0][0] - $tracer_N*$rel_mass_tracer_isotope;
      $sum_rel_masses->[$r][1] = $sum_rel_masses->[$r][1] + $merged_purity_wghts->[0][1] - $tracer_n*$rel_mass_tracer_isotope;
      #########################################################################
   }

   my $tracer_precursor = $self->{precursor}{tracer};
   my $tracer_fragment  = $self->{fragment}{tracer};
   for( my $t = 0; $t < @$merged_elems; $t++ )
   {
      if( $sum_rel_masses->[$t][0] > $rel_mass_tracer_isotope*$tracer_precursor ||
          $sum_rel_masses->[$t][1] > $rel_mass_tracer_isotope*$tracer_fragment )
      {
         splice @$compound_probs, $t, 1;
         splice @$merged_elems,   $t, 1;
         splice @$sum_rel_masses, $t, 1;
         $t--;
      }
   }

   # print "iso_names: @$iso_names\n";
}
################################################################################


################################################################################
################################################################################
sub _compute_combined_purity_Nn_pairs
{
   my $self     = shift;
   my $tracer_N = shift;
   my $tracer_n = shift;

   my $arr_Nn_pairs;
   my $arr_iso_names;
   my $arr_iso_rel_masses;

   my $tracer_elem = $self->{tracer_elem};
   my $nat_isotopes = $self->{nat_isotopes};

   my $lowest_isotope = $nat_isotopes->get_lowest_isotope($tracer_elem);

   my $atom_loss = $tracer_N - $tracer_n;

   foreach my $iso (sort $nat_isotopes->get_all_isotopes_for_element($tracer_elem))
   {
      # next if $iso eq $lowest_isotope;
      my $rel_mass = $nat_isotopes->get_massnumber_for_isotope_relative_to_lowest_isotope($iso);

      my $num_entries = 0;
      my $pairs = [];
      for( my $n = 0; $n <= $tracer_n; $n++ )
      {
         for( my $N = 0; $N <= $tracer_N; $N++ )
         {
            if( $n > $N )
            {
            }
            elsif( $atom_loss + $n < $N )
            {
            }
            else
            {
               $pairs->[$num_entries][0] = $N;
               $pairs->[$num_entries][1] = $n;
               $pairs->[$num_entries][2] = $rel_mass*$N;
               $pairs->[$num_entries][3] = $rel_mass*$n;
               $num_entries++;
            }
         }
      }
      # print "Nn pairs for purity '$iso' (tracer_N=$tracer_N tracer_n=$tracer_n):\n";
      # _print_Nn_pairs($pairs);
      push @$arr_Nn_pairs, $pairs;
      push @$arr_iso_names, 'purity_tracer_' . $iso;
      push @$arr_iso_rel_masses, $rel_mass;
   }

   my $elem_name = 'purity_tracer_' . $tracer_elem;
   my ( $merged_purity, $merged_purity_wghts ) = _merge_iso_pairs_purity($self, $elem_name, $arr_Nn_pairs, $tracer_N, $tracer_n);
   # _print_merged( $merged_purity );

   return $merged_purity, $merged_purity_wghts;
}
################################################################################


################################################################################
################################################################################
sub _print_compound
{
   my $whoami = _whoami();
   my $cmpd   = shift;
   my $iso_nm = shift;

   for( my $i = 0; $i < @$cmpd; $i++ )
   {
      print " $iso_nm->[$i]_N=$cmpd->[$i][0] $iso_nm->[$i]_n=$cmpd->[$i][1]";
   }
   print "\n";
}
################################################################################


################################################################################
################################################################################
sub print_isotope_combinations
{
   my $whoami = _whoami();
   my $self = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);

   my $merged_elems   = $self->{merged_elems};
   my $compound_probs = $self->{compound_probabilities};
   my $sum_rel_masses = $self->{sum_rel_masses};
   my $iso_names      = $self->{merged_iso_names};

   for( my $r = 0; $r < @$merged_elems; $r++ )
   {
      print "$r:";
      my $num_entries = @{$merged_elems->[$r]};
      for( my $e = 0; $e < $num_entries; $e++ )
      {
         my $iso_name = $iso_names->[$e];
         if( $iso_name eq 'tracer' )
         {
            $iso_name .= "($self->{tracer_iso})";
         }
         print " N($iso_name)=$merged_elems->[$r][$e][0] n($iso_name)=$merged_elems->[$r][$e][1]"
      }
      print " M=$sum_rel_masses->[$r][0] m=$sum_rel_masses->[$r][1] prob=$compound_probs->[$r]\n";
   }
}
################################################################################


################################################################################
################################################################################
sub _n_over_k
{
   my $n = shift;
   my $k = shift;

   my $key = $n . '_' . $k;

   return $g_cache_n_over_k->{$key} if $g_cache_n_over_k->{$key};

   my $num_cmbs = 1;
   my $n_help = $n;
   for( my $i = 0; $i < $k; $i++ )
   {
      $num_cmbs *= $n_help;
      $n_help--;
      $num_cmbs /= ($i+1);
   }

   $g_cache_n_over_k->{$key} = $num_cmbs;

   return $num_cmbs;
}
################################################################################


################################################################################
################################################################################
sub _sum_up_masses
{
   my $whoami = _whoami();
   my $self = shift;
   my $arr_sum_rel_masses = [];

   my $merged_elems = $self->{merged_elems};
   my $rel_masses   = $self->{merged_iso_rel_masses};

   for( my $r = 0; $r < @$merged_elems; $r++ )
   {
      my $num_entries = @{$merged_elems->[$r]};
      $arr_sum_rel_masses->[$r][0] = 0;
      $arr_sum_rel_masses->[$r][1] = 0;

      for( my $e = 0; $e < $num_entries; $e++ )
      {
         my $rel_mass = $rel_masses->[$e];
         # print "r=$r, e=$e @{$rel_masses}\n";
         $arr_sum_rel_masses->[$r][0] += $merged_elems->[$r][$e][0]*$rel_mass;
         $arr_sum_rel_masses->[$r][1] += $merged_elems->[$r][$e][1]*$rel_mass;
      }
      # print "$whoami: $r: sum_mass_N=$arr_sum_rel_masses->[$r][0], sum_mass_n=$arr_sum_rel_masses->[$r][1]\n";
   }

   return $arr_sum_rel_masses;
}
################################################################################


################################################################################
################################################################################
sub _merge_elems
{
   my $whoami         = _whoami();
   my $self           = shift;
   my $m_isos         = shift;
   my $m_isos_wght    = shift;
   my $elem_names     = shift;
   my $res            = [];
   my $res_weight     = [];
   my $sum_nat_ab_Nn  = [];
   my $num_entries    = 0;
   my $filtered_total = 0;
   my $elem_idx       = 0;

   # print "INFO: $whoami: entered\n";

   my $tracer_N = $self->{precursor}{tracer};
   my $tracer_n = $self->{fragment}{tracer};
   my $rel_mass_tracer_isotope = $self->{rel_mass_tracer_isotope};

   # print "tracer_N=$tracer_N tracer_n=$tracer_n rel_mass_tracer_isotope=$rel_mass_tracer_isotope\n";

   ############################################################################
   # fill result array with pairs of first isotope
   ############################################################################
   for( my $x = 0; $x < @{$m_isos->[0]}; $x++ )
   {
      $sum_nat_ab_Nn->[$num_entries][0] = 0;
      $sum_nat_ab_Nn->[$num_entries][1] = 0;

      # print "INFO: $whoami: x=$x/",scalar(@{$m_isos->[0]}),"\n";
      for( my $y = 0; $y < @{$m_isos->[0][$x]}; $y++ )
      {
         $res->[$num_entries][$y][0] = $m_isos->[0][$x][$y][0];
         $res->[$num_entries][$y][1] = $m_isos->[0][$x][$y][1];

         if( $self->{nat_ab_tracer} && ($elem_names->[$elem_idx] =~ /^nat_ab_tracer/ || $elem_names->[$elem_idx] =~ /^tracer/) )
         {
            $sum_nat_ab_Nn->[$num_entries][0] += $m_isos->[0][$x][$y][0];
            $sum_nat_ab_Nn->[$num_entries][1] += $m_isos->[0][$x][$y][1];
         }
      }
      $res_weight->[$num_entries][0] = $m_isos_wght->[0][$num_entries][0];
      $res_weight->[$num_entries][1] = $m_isos_wght->[0][$num_entries][1];
      
      $num_entries++;
   }
   $elem_idx++;
   ############################################################################

   ############################################################################
   ############################################################################
   for( my $e = 1; $e < @$m_isos; $e++ )
   {
      # print "INFO: $whoami: e=$e/",scalar(@{$m_isos}),"\n";
      #########################################################################
      # copy result arrays
      #########################################################################
      my $num_res = @$res;
      for( my $p = 1; $p < @{$m_isos->[$e]}; $p++ )
      {
         # print "INFO: $whoami:    p=$p/",scalar(@{$m_isos->[$e]}),"\n";
         for( my $r = 0; $r < $num_res; $r++ )
         {
            # if( $r%100000 == 0 )
            # {
            #    my $memGB = _get_used_memory();
            #    print "INFO: $whoami:        e=$e/",scalar(@{$m_isos})," p=$p/",scalar(@{$m_isos->[$e]})," r=$r/$num_res, mem=${memGB}GB\n";
            # }
            for( my $n = 0; $n < @{$res->[$r]}; $n++ )
            {
               my $index = $r + $p*$num_res;
               # print "$whoami: i=$i, p=$p, r=$r, n=$n, index=$index, num_res=$num_res\n";
               $res->[$r + $p*$num_res][$n][0] = $res->[$r][$n][0];
               $res->[$r + $p*$num_res][$n][1] = $res->[$r][$n][1];
            }
            $res_weight->[$r + $p*$num_res][0] = $res_weight->[$r][0];
            $res_weight->[$r + $p*$num_res][1] = $res_weight->[$r][1];
            $sum_nat_ab_Nn->[$r + $p*$num_res][0] =  $sum_nat_ab_Nn->[$r][0];
            $sum_nat_ab_Nn->[$r + $p*$num_res][1] =  $sum_nat_ab_Nn->[$r][1];
         }
      }

      #########################################################################
      # append new pair to result array
      #########################################################################
      my $r = 0;
      for( my $p = 0; $p < @{$m_isos->[$e]}; $p++ )
      {
         for( my $x = 0; $x < $num_res; $x++ )
         {
            # if( $x%100000 == 0 )
            # {
            #    my $memGB = _get_used_memory();
            #    print "INFO: $whoami:        e=$e/",scalar(@{$m_isos})," p=$p/",scalar(@{$m_isos->[$e]})," x=$x/$num_res, mem=${memGB}GB\n";
            # }

            my $nums = @{$res->[$r]};
            for( my $y = 0; $y < @{$m_isos->[$e][$p]}; $y++)
            {
               $res->[$r][$nums + $y][0] = $m_isos->[$e][$p][$y][0];
               $res->[$r][$nums + $y][1] = $m_isos->[$e][$p][$y][1];
               if( $self->{nat_ab_tracer} && ($elem_names->[$elem_idx] =~ /^nat_ab_tracer/ || $elem_names->[$elem_idx] =~ /^tracer/) )
               {
                  $sum_nat_ab_Nn->[$r][0] += $m_isos->[$e][$p][$y][0];
                  $sum_nat_ab_Nn->[$r][1] += $m_isos->[$e][$p][$y][1];
               }
            }
            $res_weight->[$r][0] += $m_isos_wght->[$e][$p][0];
            $res_weight->[$r][1] += $m_isos_wght->[$e][$p][1];
            $r++;
         }
      }
      #########################################################################


      #########################################################################
      # get rid of isotope combinations which are irrelevant
      #########################################################################
      my $filtered_this_time = 0;
      for( my $t = 0; $t < @$res; $t++ )
      {
         if( $self->{tracer_purity} )
         {
            if( ($res_weight->[$t][0] - $res->[$t][-1][0]*$rel_mass_tracer_isotope) > $tracer_N*$rel_mass_tracer_isotope ||
                ($res_weight->[$t][1] - $res->[$t][-1][1]*$rel_mass_tracer_isotope) > $tracer_n*$rel_mass_tracer_isotope )
            {
               splice @$res,           $t, 1;
               splice @$res_weight,    $t, 1;
               splice @$sum_nat_ab_Nn, $t, 1;
               $t--;
               $filtered_this_time++;
               $filtered_total++;
            }
            if( $self->{nat_ab_tracer} )
            {
               if( $sum_nat_ab_Nn->[$t][0] > $tracer_N ||
                   $sum_nat_ab_Nn->[$t][1] > $tracer_n )
               {
                  splice @$res,           $t, 1;
                  splice @$res_weight,    $t, 1;
                  splice @$sum_nat_ab_Nn, $t, 1;
                  $t--;
                  $filtered_this_time++;
                  $filtered_total++;
               }
            }
         }
         else
         {
            if( $res_weight->[$t][0] > $tracer_N*$rel_mass_tracer_isotope ||
                $res_weight->[$t][1] > $tracer_n*$rel_mass_tracer_isotope )
            {
               splice @$res,        $t, 1;
               splice @$res_weight, $t, 1;
               splice @$sum_nat_ab_Nn, $t, 1;
               $t--;
               $filtered_this_time++;
               $filtered_total++;
            }
            elsif( $self->{nat_ab_tracer} )
            {
               if( $sum_nat_ab_Nn->[$t][0] > $tracer_N ||
                   $sum_nat_ab_Nn->[$t][1] > $tracer_n )
               {
                  print "t=$t sum_nat_ab_Nn->[$t][0]=$sum_nat_ab_Nn->[$t][0] sum_nat_ab_Nn->[$t][1]=$sum_nat_ab_Nn->[$t][1] tracer_N=$tracer_N tracer_n=$tracer_n\n";
                  splice @$res,           $t, 1;
                  splice @$res_weight,    $t, 1;
                  splice @$sum_nat_ab_Nn, $t, 1;
                  $t--;
                  $filtered_this_time++;
                  $filtered_total++;
               }
            }
         }
      }
      print "INFO: $whoami: filtered_total=$filtered_total filtered_this_time=$filtered_this_time\n";
      #########################################################################
      $elem_idx++;
   }
   ############################################################################

   # print "INFO: $whoami: leaving\n";

   return $res;
}
################################################################################


################################################################################
################################################################################
sub _merge_iso_pairs
{
   my $whoami       = _whoami();
   my $self         = shift;
   my $elem         = shift;
   my $arr_Nn_pairs = shift;
   my $max_N        = shift;
   my $max_n        = shift;
   my $atom_loss    = $max_N - $max_n;
   my $res = [];
   my $res_weight = [];
   my $ret = [];
   my $ret_weight = [];
   my $num_entries = 0;
   my $num_isos    = 0;

   print "INFO: $whoami: entered: $elem, max_N=$max_N, max_n=$max_n\n";
   my $tracer_N = $self->{precursor}{tracer};
   my $tracer_n = $self->{fragment}{tracer};
   my $rel_mass_tracer_isotope = $self->{rel_mass_tracer_isotope};

   ############################################################################
   # fill result array with pairs of first isotope
   ############################################################################
   for( my $p = 0; $p < @{$arr_Nn_pairs->[0]}; $p++ )
   {
      $res->[$num_entries][$num_isos][0] = $arr_Nn_pairs->[0][$p][0];
      $res->[$num_entries][$num_isos][1] = $arr_Nn_pairs->[0][$p][1];
      $res_weight->[$num_entries][0] = $arr_Nn_pairs->[0][$p][2];
      $res_weight->[$num_entries][1] = $arr_Nn_pairs->[0][$p][3];
      $num_entries++;
   }
   $num_isos++;
   ############################################################################

   ############################################################################
   ############################################################################
   for( my $i = 1; $i < @$arr_Nn_pairs; $i++ )
   {
      #########################################################################
      # copy result arrays
      #########################################################################
      my $num_res = @$res;
      for( my $p = 1; $p < @{$arr_Nn_pairs->[$i]}; $p++ )
      {
         for( my $r = 0; $r < $num_res; $r++ )
         {
            for( my $n = 0; $n < @{$res->[$r]}; $n++ )
            {
               my $index = $r + $p*$num_res;
               # print "$whoami: i=$i, p=$p, r=$r, n=$n, index=$index, num_res=$num_res\n";
               $res->[$r + $p*$num_res][$n][0] = $res->[$r][$n][0];
               $res->[$r + $p*$num_res][$n][1] = $res->[$r][$n][1];
               $res_weight->[$r + $p*$num_res][0] = $res_weight->[$r][0];
               $res_weight->[$r + $p*$num_res][1] = $res_weight->[$r][1];
            }
         }
      }
      #########################################################################

      #########################################################################
      # append new pair to result array
      #########################################################################
      my $r = 0;
      for( my $p = 0; $p < @{$arr_Nn_pairs->[$i]}; $p++ )
      {
         for( my $x = 0; $x < $num_res; $x++ )
         {
            my $nums = @{$res->[$r]};
            $res->[$r][$nums][0]  = $arr_Nn_pairs->[$i][$p][0];
            $res->[$r][$nums][1]  = $arr_Nn_pairs->[$i][$p][1];
            $res_weight->[$r][0] += $arr_Nn_pairs->[$i][$p][2];
            $res_weight->[$r][1] += $arr_Nn_pairs->[$i][$p][3];
            $r++;
         }
      }
      #########################################################################
   }
   ############################################################################

   ############################################################################
   # filter new set of merged pairs
   ############################################################################
   my $num_excluded = 0;
   for( my $r = 0; $r < @$res; $r++ )
   {
      my $sum_N = 0;
      my $sum_n = 0;
      for( my $n = 0 ; $n < @{$res->[$r]}; $n++ )
      {
         $sum_N += $res->[$r][$n][0];
         $sum_n += $res->[$r][$n][1];
      }

      if( $res_weight->[$r][0] <= $tracer_N*$rel_mass_tracer_isotope ||
          $res_weight->[$r][1] <= $tracer_n*$rel_mass_tracer_isotope )
      {
         if( $sum_N <= $max_N && $sum_n <= $max_n && $atom_loss + $sum_n >= $sum_N )
         {
            push @$ret, $res->[$r];
            push @$ret_weight, $res_weight->[$r];
         }
         else
         {
            $num_excluded++;
         }
      }
      else
      {
         $num_excluded++;
      }
   }
   # print "INFO: $whoami: num_excluded=$num_excluded\n";
   ############################################################################

   # _print_merged($ret);

   print "INFO: $whoami: leaving\n";

   return $ret, $ret_weight;
}
################################################################################


################################################################################
################################################################################
sub _merge_iso_pairs_purity
{
   my $whoami       = _whoami();
   my $self         = shift;
   my $elem         = shift;
   my $arr_Nn_pairs = shift;
   my $max_N        = shift;
   my $max_n        = shift;
   my $atom_loss    = $max_N - $max_n;
   my $res = [];
   my $res_weight = [];
   my $ret = [];
   my $ret_weight = [];
   my $num_entries = 0;
   my $num_isos    = 0;

   # print "INFO: $whoami: entered: $elem, max_N=$max_N, max_n=$max_n\n";
   my $rel_mass_tracer_isotope = $self->{rel_mass_tracer_isotope};

   ############################################################################
   # fill result array with pairs of first isotope
   ############################################################################
   for( my $p = 0; $p < @{$arr_Nn_pairs->[0]}; $p++ )
   {
      $res->[$num_entries][$num_isos][0] = $arr_Nn_pairs->[0][$p][0];
      $res->[$num_entries][$num_isos][1] = $arr_Nn_pairs->[0][$p][1];
      $res_weight->[$num_entries][0] = $arr_Nn_pairs->[0][$p][2];
      $res_weight->[$num_entries][1] = $arr_Nn_pairs->[0][$p][3];
      $num_entries++;
   }
   $num_isos++;
   ############################################################################

   ############################################################################
   ############################################################################
   for( my $i = 1; $i < @$arr_Nn_pairs; $i++ )
   {
      #########################################################################
      # copy result arrays
      #########################################################################
      my $num_res = @$res;
      for( my $p = 1; $p < @{$arr_Nn_pairs->[$i]}; $p++ )
      {
         for( my $r = 0; $r < $num_res; $r++ )
         {
            for( my $n = 0; $n < @{$res->[$r]}; $n++ )
            {
               my $index = $r + $p*$num_res;
               # print "$whoami: i=$i, p=$p, r=$r, n=$n, index=$index, num_res=$num_res\n";
               $res->[$r + $p*$num_res][$n][0] = $res->[$r][$n][0];
               $res->[$r + $p*$num_res][$n][1] = $res->[$r][$n][1];
               $res_weight->[$r + $p*$num_res][0] = $res_weight->[$r][0];
               $res_weight->[$r + $p*$num_res][1] = $res_weight->[$r][1];
            }
         }
      }
      #########################################################################

      #########################################################################
      # append new pair to result array
      #########################################################################
      my $r = 0;
      for( my $p = 0; $p < @{$arr_Nn_pairs->[$i]}; $p++ )
      {
         for( my $x = 0; $x < $num_res; $x++ )
         {
            my $nums = @{$res->[$r]};
            $res->[$r][$nums][0]  = $arr_Nn_pairs->[$i][$p][0];
            $res->[$r][$nums][1]  = $arr_Nn_pairs->[$i][$p][1];
            $res_weight->[$r][0] += $arr_Nn_pairs->[$i][$p][2];
            $res_weight->[$r][1] += $arr_Nn_pairs->[$i][$p][3];
            $r++;
         }
      }
      #########################################################################
   }
   ############################################################################

   ############################################################################
   # filter new set of merged pairs
   ############################################################################
   my $num_excluded = 0;
   for( my $r = 0; $r < @$res; $r++ )
   {
      my $sum_N = 0;
      my $sum_n = 0;
      for( my $n = 0 ; $n < @{$res->[$r]}; $n++ )
      {
         $sum_N += $res->[$r][$n][0];
         $sum_n += $res->[$r][$n][1];
      }

      # print "r=$r: sum_N=$sum_N sum_n=$sum_n max_N=$max_N max_n=$max_n atom_loss=$atom_loss\n";
      if( $res_weight->[$r][0] <= $max_N*$rel_mass_tracer_isotope ||
          $res_weight->[$r][1] <= $max_n*$rel_mass_tracer_isotope )
      {
         # if( $sum_N <= $max_N && $sum_n <= $max_n && $atom_loss + $sum_n >= $sum_N )
         if( $sum_N == $max_N && $sum_n == $max_n && $atom_loss + $sum_n == $sum_N )
         {
            push @$ret, $res->[$r];
            push @$ret_weight, $res_weight->[$r];
         }
         else
         {
            $num_excluded++;
         }
      }
      else
      {
         $num_excluded++;
      }
   }
   print "INFO: $whoami: num_excluded=$num_excluded\n";
   ############################################################################

   # _print_merged($ret);

   # print "INFO: $whoami: leaving\n";

   return $ret, $ret_weight;
}
################################################################################


################################################################################
################################################################################
sub _print_merged
{
   my $whoami = _whoami();
   my $merged = shift;

   print "$whoami: entered.\n";
   for( my $e = 0; $e < @$merged; $e++ )
   {
      print "$whoami: $e:";
      for( my $i = 0; $i < @{$merged->[$e]}; $i++ )
      {
         print " (N=$merged->[$e][$i][0], n=$merged->[$e][$i][1])";
      }
      print "\n";
   }
   print "$whoami: leaving.\n";
}
################################################################################


################################################################################
################################################################################
sub print_tracer_Nn_pairs
{
   my $whoami    = _whoami();
   my $self      = shift;
   my $num_pairs = 0;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);

   _print_Nn_pairs($self->{tracer_Nn_pairs});
}
################################################################################


################################################################################
################################################################################
sub _print_Nn_pairs
{
   my $whoami = _whoami();
   my $Nn_pairs = shift;

   for( my $p = 0; $p < @$Nn_pairs; $p++ )
   {
      print "$p: N=$Nn_pairs->[$p][0], n=$Nn_pairs->[$p][1]\n";
   }
}
################################################################################


################################################################################
################################################################################
sub print_tracer_iso_matrix
{
   my $whoami = _whoami();
   my $self = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);

   my $mat = $self->{c_iso_matrix};

   print "    ";
   for( my $M = 0; $M < @{$mat->[0]}; $M++ )
   {
      print " M+$M";
   }
   print "\n";

   for( my $m = 0; $m < @$mat; $m++ )
   {
      print "m+$m:";
      for( my $M = 0; $M < @{$mat->[$m]}; $M++ )
      {
         print "   $mat->[$m][$M]";
      }
      print "\n";
   }
}
################################################################################


################################################################################
################################################################################
sub _compute_Nn_pairs
{
   my $whoami    = _whoami();
   my $self      = shift;
   my $elem      = shift;
   my $rel_mass  = shift;
   my $pairs     = [];
   my $skipped_p = 0;

   die "ERROR: $whoami: number of '$elem' atoms in precursor not defined" unless defined $self->{precursor}{elements}{$elem};
   die "ERROR: $whoami: number of '$elem' atoms in fragment not defined"  unless defined $self->{fragment}{elements}{$elem};

   if( $self->{precursor}{elements}{$elem} < $self->{fragment}{elements}{$elem} )
   {
      die "ERROR: $whoami: number of atoms in precursor ($self->{precursor}{elements}{$elem}) is less than in fragment ($self->{fragment}{elements}{$elem})";
   }

   my $precursor_tracer_num  = $self->{precursor}{tracer};
   my $fragment_tracer_num   = $self->{fragment}{tracer};
   my $rel_mass_tracer_isotope = $self->{rel_mass_tracer_isotope};
   my $rel_mass_tracer_isotope_precursor = $precursor_tracer_num*$rel_mass_tracer_isotope;
   my $rel_mass_tracer_isotope_fragment  = $fragment_tracer_num*$rel_mass_tracer_isotope;

   my $atom_loss = $self->{precursor}{elements}{$elem} - $self->{fragment}{elements}{$elem};

   my $num_entries = 0;
   for( my $n = 0; $n <= $self->{fragment}{elements}{$elem}; $n++ )
   {
      for( my $N = 0; $N <= $self->{precursor}{elements}{$elem}; $N++ )
      {
         if( $n > $N )
         {
         }
         elsif( $atom_loss + $n < $N )
         {
         }
         else
         {
            if( $rel_mass*$N > $rel_mass_tracer_isotope_precursor ||
                $rel_mass*$n > $rel_mass_tracer_isotope_fragment )
            {
               $skipped_p++;
            }
            else
            {
               $pairs->[$num_entries][0] = $N;
               $pairs->[$num_entries][1] = $n;
               $pairs->[$num_entries][2] = $rel_mass*$N;
               $pairs->[$num_entries][3] = $rel_mass*$n;
               $num_entries++;
            }
         }
      }
   }

   # print "INFO: $whoami: elem=$elem, rel_mass=$rel_mass: num_entries=$num_entries, skipped_p=$skipped_p\n";

   return $pairs;
}
################################################################################


################################################################################
################################################################################
sub _compute_Nn_pairs_nat_ab_tracer
{
   my $whoami    = _whoami();
   my $self      = shift;
   my $elem      = shift;
   my $rel_mass  = shift;
   my $pairs     = [];
   my $skipped_p = 0;

   # die "ERROR: $whoami: number of '$elem' atoms in precursor not defined" unless defined $self->{precursor}{elements}{$elem};
   # die "ERROR: $whoami: number of '$elem' atoms in fragment not defined"  unless defined $self->{fragment}{elements}{$elem};

   # if( $self->{precursor}{elements}{$elem} < $self->{fragment}{elements}{$elem} )
   # {
   #    die "ERROR: $whoami: number of atoms in precursor ($self->{precursor}{elements}{$elem}) is less than in fragment ($self->{fragment}{elements}{$elem})";
   # }

   my $precursor_tracer_num  = $self->{precursor}{tracer};
   my $fragment_tracer_num   = $self->{fragment}{tracer};
   my $rel_mass_tracer_isotope = $self->{rel_mass_tracer_isotope};
   my $rel_mass_tracer_isotope_precursor = $precursor_tracer_num*$rel_mass_tracer_isotope;
   my $rel_mass_tracer_isotope_fragment  = $fragment_tracer_num*$rel_mass_tracer_isotope;

   my $atom_loss = $self->{precursor}{tracer} - $self->{fragment}{tracer};

   my $num_entries = 0;
   for( my $n = 0; $n <= $self->{fragment}{tracer}; $n++ )
   {
      for( my $N = 0; $N <= $self->{precursor}{tracer}; $N++ )
      {
         if( $n > $N )
         {
         }
         elsif( $atom_loss + $n < $N )
         {
         }
         else
         {
            #if( $rel_mass*$N > $rel_mass_tracer_isotope_precursor ||
            #    $rel_mass*$n > $rel_mass_tracer_isotope_fragment )
            #{
            #   $skipped_p++;
            #}
            #else
            {
               $pairs->[$num_entries][0] = $N;
               $pairs->[$num_entries][1] = $n;
               $pairs->[$num_entries][2] = $rel_mass*$N;
               $pairs->[$num_entries][3] = $rel_mass*$n;
               $num_entries++;
            }
         }
      }
   }

   print "INFO: $whoami: elem=$elem, rel_mass=$rel_mass: num_entries=$num_entries, skipped_p=$skipped_p\n";

   return $pairs;
}
################################################################################


################################################################################
################################################################################
sub _compute_tracer_iso_matrix
{
   my $whoami       = _whoami();
   my $self         = shift;
   my $nat_isotopes = $self->{nat_isotopes};
   my $mat          = [];
   my $Nn_pairs     = [];
   my $num_pairs    = 0;

   print "INFO: $whoami: entered\n";

   # do some checks first
   die "ERROR: $whoami: number of carbon atoms in precursor not defined" unless defined $self->{precursor}{tracer};
   die "ERROR: $whoami: number of carbon atoms in fragment not defined"  unless defined $self->{fragment}{tracer};

   if( $self->{precursor}{tracer} < $self->{fragment}{tracer} )
   {
      die "ERROR: $whoami: number of tracer atoms in precursor ($self->{precursor}{tracer}) is less than in fragment ($self->{fragment}{tracer})";
   }

   my $tracer_loss = $self->{precursor}{tracer} - $self->{fragment}{tracer};
   my $rel_mass_tracer_iso= $self->{rel_mass_tracer_isotope};

   # carp "INFO: $whoami: tracer_lost=$tracer_loss fragment_C=$self->{fragment}{tracer} precursor_C=$self->{precursor}{tracer}\n";
   my $tracer_elem = $self->{tracer_elem};
   print "INFO: $whoami: tracer element '$tracer_elem'\n";
   my $tracer_iso = $self->{tracer_iso};
   print "INFO: $whoami: tracer element '$tracer_iso'\n";

   my $rel_mass_number = $nat_isotopes->get_massnumber_for_isotope_relative_to_lowest_isotope($tracer_iso);

   for( my $m = 0; $m <= $self->{fragment}{tracer}; $m++ )
   {
      for( my $M = 0; $M <= $self->{precursor}{tracer}; $M++ )
      {
         if( $m > $M )
         {
            $mat->[$m][$M] = 0;
         }
         elsif( $tracer_loss + $m < $M )
         {
            $mat->[$m][$M] = 0;
         }
         else
         {
            $mat->[$m][$M] = 1;
            $Nn_pairs->[$num_pairs][0] = $M;
            $Nn_pairs->[$num_pairs][1] = $m;
            $Nn_pairs->[$num_pairs][2] = $M*$rel_mass_tracer_iso;
            $Nn_pairs->[$num_pairs][3] = $m*$rel_mass_tracer_iso;
            $num_pairs++;
         }
      }
   }

   # check if we have the correct number of measured TOF values to be corrected
   my $num_tof_measured = @{$self->{tof_measured}};
   if( $num_tof_measured!= $num_pairs )
   {
      # die "ERROR: number of measured TOF values ($num_tof_measured) is not equal to number of non-zero elements in isotope matrix ($num_pairs)";
      carp "ERROR: number of measured TOF values ($num_tof_measured) is not equal to number of non-zero elements in isotope matrix ($num_pairs)";
   }

   print "INFO: $whoami: leaving\n";

   return $mat, $Nn_pairs;
}
################################################################################


################################################################################
################################################################################
sub get_precursor_elements
{
   my $whoami = _whoami();
   my $self = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);
   return sort keys %{$self->{precursor}};

}
################################################################################


################################################################################
################################################################################
sub get_fragment_elements
{
   my $whoami = _whoami();
   my $self = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);
   return sort keys %{$self->{fragment}};

}
################################################################################


################################################################################
################################################################################
sub get_precursor_num_entities_for_element
{
   my $whoami = _whoami();
   my $self = shift;
   my $elem = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);
   return $self->{precursor}{$elem};
}
################################################################################


################################################################################
################################################################################
sub get_fragment_num_entities_for_element
{
   my $whoami = _whoami();
   my $self = shift;
   my $elem = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);
   return $self->{fragment}{$elem};
}
################################################################################


################################################################################
# a) check if we have tof values for all M-m combinations of the tracer isoptope
# b) check if we do not have tof values which do not map on M-m combinations
################################################################################
sub _check_tof_values
{
   my $whoami = _whoami();
   my $self = shift;

   my $tracer_Nn_pairs = $self->{tracer_Nn_pairs};
   my $tof_hash        = $self->{tof_hash};

   for( my $p = 0; $p < @$tracer_Nn_pairs; $p++ )
   {
      my $M = $tracer_Nn_pairs->[$p][2];
      my $m = $tracer_Nn_pairs->[$p][3];
      my $k = $M . '_' . $m;

      unless( defined $tof_hash->{$k} )
      {
         die "ERROR: $whoami: didn't find measured value for M=$M and m=$m, check your input file";
      }
   }

   foreach my $key (sort keys %{$tof_hash})
   {
      my ($M, $m);
      if( $key =~ /^(\d+)_(\d+)$/ )
      {
         $M = $1;
         $m = $2;
      }
      else
      {
         die "ERROR: $whoami: invalid key '$key', expected: number, underscore, and number";
      }

      my $found_pair = 0;

      for( my $p = 0; $p < @$tracer_Nn_pairs; $p++ )
      {
         my $N = $tracer_Nn_pairs->[$p][2];
         my $n = $tracer_Nn_pairs->[$p][3];

         if( $M == $N && $m == $n )
         {
            $found_pair = 1;
            last;
         }
      }

      if( $found_pair == 0 )
      {
         die "ERROR: $whoami: measured value for M=$M and m=$m does not makes sense, is not relevant for this input data";
      }
   }

   # if we end up here, everything seems to be OK
}
################################################################################


################################################################################
################################################################################
sub _read_tof_values
{
   my $whoami = _whoami();
   my $fh   = shift;
   my $fn   = shift;
   my $self = shift;
   my $tof_measured;
   my $tof_hash;

   while(my $line = <$fh> )
   {
      chomp $line;
      next if $line =~ /^\s*#/;
      $line =~ s/\s+//g;

      my ($prefix, $val) = split ',', $line;

      my ($comp, $N, $n);
      if( $prefix =~ /^([_A-Za-z0-9\-]+)_M?(\d+)\.m?(\d+)$/ )
      {
         $comp = $1;
         $N    = $2;
         $n    = $3;
      }
      elsif( $prefix =~ /^([_A-Za-z0-9\-]+)_M?(\d+)$/ )
      {
         $comp = $1;
         $N    = $2;
         $n    = $2;
      }
      else
      {
         die "ERROR: $whoami: invalid prefix '$prefix' in line $. '$line'";
      }

      my $entry = {};
      $entry->{N}   = $N;
      $entry->{n}   = $n;
      $entry->{val} = $val;
    
      push @$tof_measured, $entry;

      my $key = $N . '_' . $n;
      if( exists $tof_hash->{$key} )
      {
         die "ERROR: $whoami: there was already a TOF value for (N=$N,n=$n)";
      }
      $tof_hash->{$key} = $val;
   }

   $self->{tof_measured} = $tof_measured;
   $self->{tof_hash}     = $tof_hash;
}
################################################################################


################################################################################
################################################################################
sub _read_purity_info
{
   my $whoami = _whoami();
   my $fh   = shift;
   my $fn   = shift;
   my $self = shift;

   my $nat_isotopes = $self->{nat_isotopes};

   my $line;
   while($line = <$fh>)
   {
      # get rid of newline character at the end of line
      chomp $line;
      next if $line =~ /^\s*$/;
      next if $line =~ /^\s*#/;
      last;
   }

   # trim line
   $line =~ s/^\s*//;
   $line =~ s/\s*$//;

   my ($str_iso_names, $str_purity) = split /:/, $line;

   # trim strings
   $str_iso_names =~ s/^\s*//;
   $str_iso_names =~ s/\s*$//;
   $str_purity    =~ s/^\s*//;
   $str_purity    =~ s/\s*$//;

   die "ERROR: $whoami: invalid syntax in purity file '$fn': '$line'\n" unless $str_iso_names && $str_purity;

   my @iso_names    = split /\s+/, $str_iso_names;
   my @purities     = split /\s+/, $str_purity;
   my $num_isos     = @iso_names;
   my $num_puris    = @purities;
   my $sum_purities = 0;

   die "ERROR: $whoami: invalid syntax in purity file '$fn': number of isotopes ($num_isos) and number of purities ($num_puris) differ\n" unless $num_isos == $num_puris;

   my $tracer_iso  = $self->{tracer_iso};
   my $tracer_elem = $self->{tracer_elem};

   my %all_isos_of_tracer_elem = map {$_ => 0} $nat_isotopes->get_all_isotopes_for_element($tracer_elem);

   for( my $i = 0; $i < $num_isos; $i++ )
   {
      if( $iso_names[$i] =~ /^([A-Z][a-z]?)(\d+)$/ )
      {
         my $t_elem = $1;
         if( $tracer_elem ne $t_elem )
         {
            die "ERROR: $whoami: invalid isotope '$iso_names[$i]' name in purity file '$fn': tracer element is $tracer_elem and not $t_elem\n";
         }

         die "ERROR: $whoami: isotope $iso_names[$i] found in purity file '$fn' was not found in isotope data set\n" unless  exists $all_isos_of_tracer_elem{$iso_names[$i]};

         $all_isos_of_tracer_elem{$iso_names[$i]}++;
      }
      else
      {
         die "ERROR: $whoami: invalid isotope '$iso_names[$i]' name in purity file '$fn'\n";
      }
      $self->{purity}{$iso_names[$i]} = $purities[$i];
      $sum_purities += $purities[$i];
   }

   if( $sum_purities != 1.0 )
   {
      die "ERROR: $whoami: sum of purities ($sum_purities) is not equal to 1\n";
   }

   foreach my $key (sort keys %all_isos_of_tracer_elem )
   {
      die "ERROR: $whoami: we found isoptope '$key' $all_isos_of_tracer_elem{$key} times in file '$fn'. Expected one occurrence.\n" if $all_isos_of_tracer_elem{$key} != 1;
   }
}
################################################################################


################################################################################
################################################################################
sub _process_data_sets
{
   my $whoami = _whoami();
   my $fh   = shift;
   my $fn   = shift;
   my $self = shift;
   my $found_precursor = 0;
   my $found_fragment  = 0;
   my $nat_isotopes    = $self->{nat_isotopes};

   while(my $line = <$fh>)
   {
      # get rid of newline character at the end of line
      chomp $line;

      # do not process comment lines
      next if $line =~ /^\s*#/;

      # do not process empty lines
      next if $line =~ /^\s*$/;

      $line =~ s/\s+//g;

      # extract tag of data line
      my ($tag, $data) = split ':', $line;
      die "ERROR: $whoami: Invalid line '$line' found in file '$fn' containing chemical data." unless defined $tag && defined $data;

      my ($tracer, @elems_occ) = split ',', $data;

      if( lc $tag eq 'precursor' )
      {
         $tag = 'precursor';
         $found_precursor = 1;
      }
      elsif( lc $tag eq 'production' )
      {
         $tag = 'fragment';
         $found_fragment = 1;
      }
      else
      {
         die "ERROR: $whoami: Invalid tag '$tag' found in file '$fn'. Allowed tags: 'Precursor' and 'Product Ion'";
      }

      #########################################################################
      #########################################################################
      if( $tracer =~ /^(\d+)([A-Z][a-z]?)(\d+)$/ )
      {
         my $tracer_iso = $2 . $1;
         if( defined $self->{tracer_elem} && $self->{tracer_elem} ne $2 )
         {
            die "ERROR: $whoami: tracer element of precursor and of fragment are not identical ($self->{tracer_elem}/$2)\n";
         }
         if( defined $self->{tracer_iso} && $self->{tracer_iso} ne $tracer_iso )
         {
            die "ERROR: $whoami: tracer isotope of precursor and of fragment are not identical ($self->{tracer_iso}/$tracer_iso)\n";
         }
         $self->{tracer_iso}   = $tracer_iso;
         $self->{tracer_elem}  = $2;
         $self->{$tag}{tracer} = $3;

         my $num_C_isotopes = $nat_isotopes->get_number_of_isotopes($self->{tracer_elem});
         die "ERROR: $whoami: Invalid number of isotopes ($num_C_isotopes) for tracer/tracer element. Expected more than one isotope\n" unless $num_C_isotopes > 1;

         my $rel_mass_tracer_isotope = $nat_isotopes->get_massnumber_for_isotope_relative_to_lowest_isotope($tracer_iso);
         $self->{rel_mass_tracer_isotope} = $rel_mass_tracer_isotope;
         print "INFO: $whoami: rel_mass_tracer_isotope=$rel_mass_tracer_isotope\n";
      }
      else
      {
         die "ERROR: $whoami: Invalid chemical element for '$tag' found in file '$fn': '$tracer'";
      }
      #########################################################################

      for( my $i = 0; $i < @elems_occ; $i++ )
      {
         if( $elems_occ[$i] =~ /^([A-Z][a-z]?)(\d+)$/ )
         {
            $self->{$tag}{elements}{$1} = $2;
         }
         else
         {
            die "ERROR: $whoami: Invalid chemical element for '$tag' found in file '$fn': '$elems_occ[$i]'";
         }
      }
   }

   # check if each element has a counterpart
   my %seen;
   # foreach my $elem (keys %{$self->{precursor}{elements}}, keys %{$self->{fragment}{elements}} )
   foreach my $elem (keys(%{$self->{precursor}{elements}}), keys(%{$self->{fragment}{elements}}) )
   {
      $seen{$elem}++;
   }
   foreach my $elem (keys %seen)
   {
      die "ERROR: $whoami: element '$elem' either not found in 'precursor' or not in 'fragment' in file '$fn'" unless $seen{$elem} == 2;
   }

   die "ERROR: $whoami: No chemical data for precursor provided in file '$fn'. Use tag 'precursor'." unless $found_precursor;
   die "ERROR: $whoami: No chemical data for fragment provided in file '$fn'. Use tag 'fragment'"    unless $found_fragment;

   # check if it is an isotopologue
   # which means that precursor and fragment have the same chemical structure
   my $is_isotopologue = 1;
   foreach my $elem ( keys(%{$self->{precursor}{elements}}) )
   {
      if( $self->{precursor}{elements}{$elem} != $self->{fragment}{elements}{$elem} )
      {
         $is_isotopologue = 0;
         last;
      }
   }
   $is_isotopologue = 0 if $self->{precursor}{tracer} != $self->{fragment}{tracer};
   print "INFO: $whoami: is_isotopologue: $is_isotopologue\n";

   # move all elements that have only one isotope to 'excluded'-tag
   my $exclude_elems;
   foreach my $elem ( keys(%{$self->{precursor}{elements}}) )
   {
      print "DEBUG: $whoami: exluding element: $elem\n";
      if( $nat_isotopes->get_number_of_isotopes($elem) == 1 )
      {
         $self->{excluded_precursor}{elements}{$elem} = $self->{precursor}{elements}{$elem};
         $self->{excluded_fragment}{elements}{$elem}  = $self->{fragment}{elements}{$elem};
         push @$exclude_elems, $elem;
      }
   }

   foreach my $elem_to_exclude (@$exclude_elems )
   {
      delete $self->{precursor}{elements}{$elem_to_exclude};
      delete $self->{fragment}{elements}{$elem_to_exclude};
   }
   
   $self->{isotopologue} = $is_isotopologue;
}
################################################################################


################################################################################
################################################################################
sub get_rel_mass_tracer_isotope
{
   my $whoami = _whoami();
   my $self = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);
   return $self->{rel_mass_tracer_isotope};
}
################################################################################


################################################################################
################################################################################
sub _get_used_memory
{
   my $whoami = _whoami();
   my $ram_gb = 0;

   if( $^O =~ /linux/i )
   {
      # we only know how to determine amount of used RAM in linux environment
      my $filename = "/proc/$$/statm";
      open my $fh , $filename or die "FATAL ERROR: Unable to open stat file '$filename': $!\n";
      my @stat = split /\s+/ , <$fh>;
      close $fh;

      my $page_size = qx 'getconf PAGESIZE';

      # print "INFO: $whoami: stat: @stat\n";

      # in GByte
      # return int(1000*$stat[0]*$page_size/1024/1024/1024)/1000; # virtual memory?
      $ram_gb = int(1000*$stat[1]*$page_size/1024/1024/1024)/1000; 
   }
   else
   {
      print "INFO: we are not in a Linux environment -> can't determine memory usage\n";
   }

   return $ram_gb;
}
################################################################################

################################################################################
################################################################################
sub _print_matrix
{
   my $whoami = _whoami();
   my $mat    = shift;
   my $msg    = shift;

   print $msg;

   for( my $r = 0; $r < @$mat; $r++ )
   {
      for( my $c = 0; $c < @{$mat->[$r]}; $c++ )
      {
         print " $mat->[$r][$c]";
      }
      print "\n";
   }
}
################################################################################


###############################################################################
###############################################################################
sub _whoami
{
   ( caller(1) )[3]
}
###############################################################################

1;
