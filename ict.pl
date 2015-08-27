#! /usr/bin/perl
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

use Getopt::Std;
use File::Basename;
use lib dirname($0);
use Bio::IsotopeCorrection::NaturalIsotopes;
use Bio::IsotopeCorrection::Chemicals;
use constant VERSION => 0.04;
use constant AUTHOR  => 'Christian Jungreuthmayer';
use constant YEAR    => '2014 - 2015';

select( (select(STDOUT), $| = 1)[0] );
select( (select(STDERR), $| = 1)[0] );

# MAX_DIFF is used to determine when
# the difference between computed and
# expected result is too large and
# a warning message is printed
use constant MAX_DIFF => 1.0e-0;

use vars qw($opt_c $opt_i $opt_m $opt_e $opt_o $opt_h $opt_k $opt_n $opt_p);

my $g_filename_chemicals;
my $g_filename_nat_isotopes;
my $g_filename_measured_tof;
my $g_filename_expected_tof;
my $g_filename_corrected_output = '';
my $g_filename_gen_chemicals    = 'tmp_chem';
my $g_filename_gen_measured_tof = 'tmp_tof';
my $g_filename_purity_info      = '';
my $g_keep_tmp_files = 0;
my $g_result = [];
my $g_fho;
my $g_num_corr_proc = 0;
my $g_nat_ab_tracer = 0;

read_arguments();

my $natIsos = Bio::IsotopeCorrection::NaturalIsotopes->new({filename => $g_filename_nat_isotopes});

my $chem_data = parse_chem_file($g_filename_chemicals);
my ($tof_data, $comp_key_ind, $comp_ind)  = parse_tof_file($g_filename_measured_tof);

open_result_file( $g_filename_corrected_output ) if $g_filename_corrected_output;

# perform correction for each compound provided in input file
for( my $c = 0; $c < @$comp_key_ind; $c++ )
{
   my $comp_key = $comp_key_ind->[$c];
   my $comp     = $comp_ind->[$c];
   my $num_tof_experiments = @{$tof_data->{$comp_key}[0]{vals}};
   print "Do correction for compound $c '$comp_key'. Number of tof experiments found: $num_tof_experiments.\n";

   my $start_index_result = @$g_result;
   create_prefix($comp_key, $tof_data, $g_result) if $g_filename_corrected_output;

   # loop over all tof experiments provided for given compound
   for( my $e = 0; $e < $num_tof_experiments; $e++ )
   {
      $g_num_corr_proc++;
      my $tof_filename = $g_filename_gen_measured_tof . '_' . $comp_key . '_' . $e . '.txt';
      my $isotopologe = create_tof_file( $tof_filename, $comp_key, $tof_data, $e );

      # create information about chemical structure of compound
      # before and after fragmentation
      my $chem_filename = $g_filename_gen_chemicals . '_' . $comp_key . '_' . $e .  '.txt';
      create_chem_file($chem_filename, $comp, $chem_data, $isotopologe);

      my $chemicals = Bio::IsotopeCorrection::Chemicals->new({filename_chem => $chem_filename,
                                                              filename_tof  => $tof_filename,
                                                              nat_isotopes  => $natIsos,
                                                              nat_ab_tracer => $g_nat_ab_tracer,
                                                              tracer_purity => $g_filename_purity_info});

      # print internal information - mainly for analyzing and debugging
      # $chemicals->print_tracer_iso_matrix();
      # $chemicals->print_target_Nn_pairs();
      # $chemicals->print_isotope_combinations();

      $chemicals->do_correction();
      $chemicals->print_result();
      $chemicals->print_enrichment() if $chemicals->is_isotopologue();
      my $corrected_tof = $chemicals->get_correct_areas();

      add_corrected_tof_to_results( $g_result, $corrected_tof, $start_index_result );

      # check results against expected values provided in a separate file
      if( $g_filename_expected_tof )
      {
         check_results( $g_filename_expected_tof, $corrected_tof, $comp_key, $e );
      }

      # clean up
      unlink $chem_filename unless $g_keep_tmp_files;
      unlink $tof_filename  unless $g_keep_tmp_files;
   }

   print_result($g_filename_corrected_output, $g_result, $start_index_result) if $g_filename_corrected_output
}

close_result_file() if $g_filename_corrected_output;
print "Number of performed correction procedures: $g_num_corr_proc\n";
################################################################################


################################################################################
################################################################################
sub add_corrected_tof_to_results
{
   my $res      = shift;
   my $corr_tof = shift;
   my $start_i  = shift;

   for( my $i = 0; $i < @$corr_tof; $i++ )
   {
      push @{$res->[$start_i+$i]}, $corr_tof->[$i];
   }
}
################################################################################


################################################################################
################################################################################
sub open_result_file
{
   my $whoami   = whoami();
   my $filename = shift;

   open $g_fho, ">$filename" or die "ERROR: $whoami: couldn't open file '$filename' for writing: $!\n";

   select( (select($g_fho), $| = 1)[0] );
}
################################################################################


################################################################################
################################################################################
sub close_result_file
{
   close $g_fho;
}
################################################################################


################################################################################
################################################################################
sub print_result
{
   my $whoami   = whoami();
   my $filename = shift;
   my $result   = shift;
   my $start_i  = shift;


   for( my $i = $start_i; $i < @$result; $i++ )
   {
      print $g_fho join(', ', @{$result->[$i]}), "\n";
   }

   print "INFO: results were written to file '$filename'\n";
}
################################################################################


################################################################################
################################################################################
sub create_prefix
{
   my $whoami = whoami();
   my $comp_key    = shift;
   my $tof_info    = shift;
   my $res_prefix  = shift;

   my $start_index = @$res_prefix;

   for( my $i = 0; $i < @{$tof_info->{$comp_key}}; $i++ )
   {
      my $entry = $tof_info->{$comp_key}[$i];
      my $N = $entry->{N};
      my $n = $entry->{n};
      if( $entry->{iso} )
      {
         $res_prefix->[$start_index + $i][0] = $comp_key . '_M' . $N;
      }
      else
      {
         $res_prefix->[$start_index + $i][0] = $comp_key . '_M' . $N . '.' . $n;
      }
   }

   return $res_prefix;
}
################################################################################


################################################################################
################################################################################
sub create_tof_file
{
   my $whoami = whoami();
   my $filename    = shift;
   my $comp_key    = shift;
   my $tof_info    = shift;
   my $exp_index   = shift;
   my $isotopologe = 0;

   #############################################################################
   # open file
   #############################################################################
   open my $fh, ">$filename" or die "ERROR: $whoami: couldn't open file '$filename' for writing: $!\n";
   #############################################################################

   for( my $i = 0; $i < @{$tof_info->{$comp_key}}; $i++ )
   {
      my $entry = $tof_info->{$comp_key}[$i];
      my $comp  = $entry->{comp};
      my $N     = $entry->{N};
      my $n     = $entry->{n};
      my $v     = $entry->{vals}[$exp_index];
      $isotopologe = 1 if $entry->{iso};
      print $fh $comp . '_M' . $N . '.' . $n . ', ' . $v . "\n";
   }

   #############################################################################
   # close file
   #############################################################################
   close $fh;
   #############################################################################

   return $isotopologe;
}
################################################################################


################################################################################
################################################################################
sub create_chem_file
{
   my $whoami      = whoami();
   my $filename    = shift;
   my $comp        = shift;
   my $chem_info   = shift;
   my $isotopologe = shift;

   #############################################################################
   # open file
   #############################################################################
   open my $fh, ">$filename" or die "ERROR: $whoami: couldn't open file '$filename' for writing: $!\n";
   #############################################################################

   die "ERROR: $whoami: couldn't find chemical data for compound '$comp' -> check input files\n" unless defined $chem_info->{compounds}{$comp};

   #############################################################################
   # print data for precursor
   #############################################################################
   print $fh "Precursor:";
   print $fh ' ' .
             $chem_info->{compounds}{$comp}{precursor}{tracer_iso}  .
             $chem_info->{compounds}{$comp}{precursor}{tracer_elem} .
             $chem_info->{compounds}{$comp}{precursor}{tracer};
   foreach my $elem (keys(%{$chem_info->{compounds}{$comp}{precursor}{elements}}))
   {
      my $num = $chem_info->{compounds}{$comp}{precursor}{elements}{$elem};
      print $fh ", $elem$num";
   }
   print $fh "\n";
   #############################################################################

   #############################################################################
   # data for fragment
   #############################################################################
   print $fh "Product Ion:";
   if( $isotopologe )
   {
      # prefix in tof data was of type 'XYZ_Md' ->
      # molecule was NOT fragmented ->
      # we use data of precursor in fragment
      print $fh ' ' .
                $chem_info->{compounds}{$comp}{precursor}{tracer_iso}  .
                $chem_info->{compounds}{$comp}{precursor}{tracer_elem} .
                $chem_info->{compounds}{$comp}{precursor}{tracer};
      foreach my $elem (keys(%{$chem_info->{compounds}{$comp}{precursor}{elements}}))
      {
         my $num = $chem_info->{compounds}{$comp}{precursor}{elements}{$elem};
         print $fh ", $elem$num";
      }
   }
   else
   {
      # prefix in tof data was of type 'XYZ_Md.d' ->
      # molecule WAS fragmented ->
      # we have to use fragment data
      print $fh ' ' .
                $chem_info->{compounds}{$comp}{fragment}{tracer_iso}  .
                $chem_info->{compounds}{$comp}{fragment}{tracer_elem} .
                $chem_info->{compounds}{$comp}{fragment}{tracer};
      foreach my $elem (keys(%{$chem_info->{compounds}{$comp}{fragment}{elements}}))
      {
         my $num = $chem_info->{compounds}{$comp}{fragment}{elements}{$elem};
         print $fh ", $elem$num";
      }
   }
   print $fh "\n";
   #############################################################################

   #############################################################################
   # close file
   #############################################################################
   close $fh;
   #############################################################################
}
################################################################################


################################################################################
################################################################################
sub parse_tof_file
{
   my $whoami = whoami();
   my $filename = shift;
   my $tof_info = {};
   my $comp_key_index = [];
   my $comp_index = [];
   my $cur_key = '';

   open my $fh, $filename or die "ERROR: $whoami: couldn't open file '$filename' for reading: $!\n";

   while(my $line = <$fh>)
   {
      # get rid of newline character at the end of line
      chomp $line;

      # do not process comment lines
      next if $line =~ /^\s*#/;

      # do not process empty lines
      next if $line =~ /^\s*$/;

      $line =~ s/\s+//g;

      while( $line =~ /,,/ )
      {
         $line =~ s/,,/,0,/;
      }
      $line =~ s/,$/,0/;

      my ($prefix, @vals) = split ',', $line;

      for( my $i = 0; $i < @vals; $i++ )
      {
         $vals[$i] = 0 if $vals[$i] eq '';
      }

      my ($comp, $num, $N, $n, $isotopologe);
      if( $prefix =~ /^([A-Za-z0-9\-]+?)(_\d+)_M(\d+)\.m?(\d+)$/ )
      {
         $comp = $1;
         $num  = $2;
         $N    = $3;
         $n    = $4;
         $isotopologe = 0;
      }
      elsif( $prefix =~ /^([A-Za-z0-9\-]+?)_M(\d+)\.m?(\d+)$/ )
      {
         $comp = $1;
         $num  = '';
         $N    = $2;
         $n    = $3;
         $isotopologe = 0;
      }
      elsif( $prefix =~ /^([A-Za-z0-9\-]+?)(_\d+)_M(\d+)$/ )
      {
         $comp = $1;
         $num  = $2;
         $N    = $3;
         $n    = $N;
         $isotopologe = 1;
      }
      elsif( $prefix =~ /^([A-Za-z0-9\-]+?)_M(\d+)$/ )
      {
         $comp = $1;
         $num  = '';
         $N    = $2;
         $n    = $N;
         $isotopologe = 1;
      }
      else
      {
         die "ERROR: $whoami: invalid prefix '$prefix' in line $. '$line'";
      }

      # print "DEBUG: prefix: $prefix\n";
      # print "DEBUG: comp='$comp', num='$num', N='$N', n='$n', isotopologe='$isotopologe'\n";

      die "ERROR: $whoami: did find measured values in line $. '$line'" unless @vals;

      my $entry = {};
      $entry->{comp} = $comp;
      $entry->{N}    = $N;
      $entry->{n}    = $n;
      $entry->{num}  = $num;
      $entry->{vals} = \@vals;
      $entry->{iso}  = $isotopologe;

      my $comp_key = $comp . $num;

      # push @{$tof_info->{$comp}}, $entry;
      push @{$tof_info->{$comp_key}}, $entry;

      # if( $cur_comp ne $comp )
      if( $cur_key ne $comp_key )
      {
         # $cur_comp = $comp;
         $cur_key = $comp_key;
         push @{$comp_key_index}, $comp_key;
         push @{$comp_index}, $comp;
      }
   }

   close $fh;

   return $tof_info, $comp_key_index, $comp_index;
}
################################################################################


################################################################################
################################################################################
sub parse_chem_file
{
   my $whoami = whoami();
   my $filename = shift;
   my $chem_info = {};
   my $found_precursor = 0;
   my $found_fragment  = 0;
   my $compound_name;

   open my $fh, $filename or die "ERROR: couldn't open file '$filename' for reading: $!\n";

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
      my ($prefix, $tag, $data) = split ':', $line;
      die "ERROR: $whoami: Invalid line '$line' found in file '$filename' containing chemical data." unless defined $prefix && defined $tag && defined $data;

      my ($tracers, @elems_occ) = split ',', $data;

      if( lc($tag) eq 'precursor' )
      {
         $tag = 'precursor';
         $found_precursor = 1;

         if( $found_fragment == 1 )
         {
            if( $prefix ne $compound_name )
            {
               die "ERROR: $whoami: compound '$prefix' not equal to '$compound_name'";
            }
         }
         else
         {
            $compound_name = $prefix;
         }
      }
      elsif( lc($tag) eq 'production' )
      {
         $tag = 'fragment';
         $found_fragment = 1;

         if( $found_precursor == 1 )
         {
            if( $prefix ne $compound_name )
            {
               die "ERROR: $whoami: compound '$prefix' not equal to '$compound_name'";
            }
         }
         else
         {
            $compound_name = $prefix;
         }
      }
      else
      {
         die "ERROR: $whoami: Invalid tag '$tag' found in file '$filename'. Allowed tags: 'Precursor' and 'Product Ion'";
      }

      if( $tracers =~ /^(\d+)([A-Z][a-z]?)(\d+)$/ )
      {
         $chem_info->{compounds}{$compound_name}{$tag}{tracer_iso}  = $1;
         $chem_info->{compounds}{$compound_name}{$tag}{tracer_elem} = $2;
         $chem_info->{compounds}{$compound_name}{$tag}{tracer}      = $3;
      }
      else
      {
         die "ERROR: $whoami: Invalid chemical element for '$tag' found in file '$filename': '$tracers'";
      }

      for( my $i = 0; $i < @elems_occ; $i++ )
      {
         if( $elems_occ[$i] =~ /^([A-Z][a-z]*)(\d+)$/ )
         {
            $chem_info->{compounds}{$compound_name}{$tag}{elements}{$1} = $2;
         }
         else
         {
            die "ERROR: $whoami: Invalid chemical element for '$tag' found in file '$filename': '$elems_occ[$i]'";
         }
      }

      if( $found_fragment == 1 && $found_precursor == 1 )
      {
         $found_fragment  = 0;
         $found_precursor = 0;
         $compound_name   = '';
      }
   }
   close $fh;

   foreach my $comp (keys %{$chem_info->{compounds}} )
   {
      # check if each element has a counterpart
      my %seen;

      foreach my $elem (keys(%{$chem_info->{compounds}{$comp}{precursor}{elements}}), keys(%{$chem_info->{compounds}{$comp}{fragment}{elements}}) )
      {
         $seen{$elem}++;
      }
      foreach my $elem (keys %seen)
      {
         unless( $seen{$elem} == 2 )
         {
            die "ERROR: $whoami: compound '$comp': element '$elem' either not found in 'precursor' or not found in 'fragment' in file '$filename'\n";
         }
      }
   }



   return $chem_info;
}
################################################################################




################################################################################
################################################################################
sub check_results
{
   my $whoami     = whoami();
   my $file       = shift;
   my $cor_tof    = shift;
   my $comp       = shift;
   my $exp_ind    = shift;
   my $line_n     = 0;
   my $found_comp = 1;

   open my $fh, $file or die "ERROR: couldn't open file '$file' for reading: $!\n";

   while( my $e_tof = <$fh> )
   {
      chomp $e_tof;

      my ($read_comp, $str_vals) = split ':', $e_tof;

      next unless $read_comp eq $comp;

      $found_comp = 1;

      my @vals = split ',', $str_vals;

      # check if we have enough expected values we can test our results against
      if( $exp_ind > scalar(@vals) - 1 )
      {
         warn "ERROR: $whoami: number values  in expected results in file '$file' for compound '$comp' is invalid\n";
         die  "                exp_ind=$exp_ind, number of values ", scalar(@vals), "\n";
      }

      # compare corrected with expected results
      if( abs($vals[$exp_ind] - $cor_tof->[$line_n]) > MAX_DIFF )
      {
         warn "WARNING: expected mass spectrum value ($vals[$exp_ind]) and computed mass spectrum value ($cor_tof->[$line_n]) differ\n";
         # die "WARNING: expected mass spectrum value ($e_tof) and computed mass spectrum value ($cor_tof->[$line_n]) differ\n";
      }
      else
      {
         warn "INFO: expected mass spectrum value ($vals[$exp_ind]) and computed mass spectrum value ($cor_tof->[$line_n]) are OK\n";
      }
      $line_n++;
   }

   close $fh;

   if( $found_comp == 0 )
   {
      warn "ERROR: $whoami: we couldn't find compound '$comp' in file '$file'\n";
   }
}
################################################################################


################################################################################
################################################################################
sub read_arguments
{
   getopts('c:o:i:m:e:knp:h');

   if( $opt_h )
   {
      usage();
   }

   if( $opt_k )
   {
      $g_keep_tmp_files = 1;
   }

   if( $opt_n )
   {
      $g_nat_ab_tracer = 1;
   }

   if( $opt_i )
   {
      $g_filename_nat_isotopes = $opt_i;
   }

   if( $opt_e )
   {
      $g_filename_expected_tof = $opt_e;
   }

   if( $opt_o )
   {
      $g_filename_corrected_output = $opt_o;
   }

   if( $opt_p )
   {
      $g_filename_purity_info = $opt_p;
   }

   if( $opt_c )
   {
      $g_filename_chemicals = $opt_c;
   }
   else
   {
      usage('ERROR: Name of input file containing chemical information about precursor and fragment not provided',-1);
   }

   if( $opt_m )
   {
      $g_filename_measured_tof = $opt_m;
   }
   else
   {
      usage('ERROR: Name of input file containing measured mass spectrum values not provided',-1);
   }
}
################################################################################


################################################################################
################################################################################
sub whoami
{  
   ( caller(1) )[3]
}
################################################################################


################################################################################
################################################################################
sub get_program_name
{
   return fileparse($0);
}
################################################################################


################################################################################
################################################################################
sub usage
{
   my $message   = shift || '';
   my $exit_code = shift || 0;

   print "$message\n" if $message;

   my $prog_name = get_program_name();

   print "$prog_name -c chem_data -m measured_tof [-o output_file -i nat_isotope_data -e expected_tof -k -n -h -p purity_file]\n";
   print "\n";
   print "Version: " . VERSION . "\n";
   print "Author:  " . AUTHOR  . "\n";
   print "Year: "    . YEAR    . "\n";
   print "\n";
   print "-c ..... name of input file containing chemical information about precursor and fragment\n";
   print "-m ..... name of input file containing measured mass spectrum values\n";
   print "-o ..... name of output file corrected mass spectrum values are written to\n";
   print "-i ..... name of input file containing chemical/physical information about natural isotopes\n";
   print "-e ..... name of input file containing expected corrected mass spectrum values (used for checking this program)\n";
   print "-k ..... keep (do not remove) temporary files\n";
   print "-p ..... name of input file containing information about purity of tracer\n";
   print "-n ..... natural abundance of tracer (of not labeled fraction)\n";
   print "-h ..... print this message\n";

   exit $exit_code;
}
################################################################################
