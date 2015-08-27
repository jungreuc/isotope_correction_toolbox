package Bio::IsotopeCorrection::NaturalIsotopes;
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

our $VERSION = '0.01';

sub new
{
   my $whoami   = _whoami();
   my $class    = shift;
   my $param    = shift;
   my $filename;
   my $self = {};

   die "ERROR: $whoami: It seems that new() was already called for this object" if ref $class;

   $filename = $param->{filename} if defined $param->{filename};

   if( defined $filename )
   {
      my $num_sets_external = 0;
      open my $fh, $filename or die "ERROR: $whoami: couldn't open file '$filename' for reading: $!";
      while( my $data_set = <$fh> )
      {
         chomp $data_set;
         next if $data_set =~ /^\s*$/;
         next if $data_set =~ /^\s*#/;
         _process_data_sets($data_set, $self);
         $num_sets_external++;
      }
      close $fh;
      warn "INFO: $whoami: Found $num_sets_external isotope data sets in file '$filename'";
   }
   else
   {
      my $num_sets_internal = 0;
      # read entire isotope data sets
      while( my $data_set = <DATA> )
      {
         chomp $data_set;
         next if $data_set =~ /^\s*$/;
         next if $data_set =~ /^\s*#/;
         _process_data_sets($data_set, $self);
         $num_sets_internal++;
      }
      warn "INFO: $whoami: Found $num_sets_internal internal isotope data sets";
   }

   $self = bless $self, $class;

   return $self;
}
################################################################################


################################################################################
################################################################################
sub get_lowest_isotope
{
   my $whoami  = _whoami();
   my $self    = shift;
   my $element = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);
   die "ERROR: $whoami: element was not provided" unless defined $element;

   die "ERROR: $whoami: element '$element' not found in data set" unless defined $self->{$element};

   return $self->{$element}{lowest_isotope};
}
################################################################################


################################################################################
################################################################################
sub get_all_elements
{
   my $whoami = _whoami();
   my $self = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);
   return sort keys %$self;
}
################################################################################


################################################################################
################################################################################
sub get_all_isotopes_for_element
{
   my $whoami  = _whoami();
   my $self    = shift;
   my $element = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);
   die "ERROR: $whoami: element was not provided" unless defined $element;

   return sort keys %{$self->{$element}{isotopes}};
}
################################################################################


################################################################################
################################################################################
sub get_massnumber_for_isotope_relative_to_lowest_isotope
{
   my $whoami = _whoami();
   my $self = shift;
   my $isotope = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);

   die "ERROR: $whoami: isotope was not defined" unless defined $isotope;

   my $element = $self->get_element_for_isotope($isotope);
   my $lowest_isotope_isotope = $self->get_lowest_isotope($element);

   my $rel_mass_number = $self->get_massnumber_for_isotope($isotope) -
                         $self->get_massnumber_for_isotope($lowest_isotope_isotope);

   return $rel_mass_number;
}
################################################################################


################################################################################
################################################################################
sub get_massnumber_for_isotope
{
   my $whoami  = _whoami();
   my $self    = shift;
   my $isotope = shift;
   my $elem;
   my $massnumber;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);

   if( $isotope =~ /^([A-Z][a-z]*)(\d+)/ )
   {
      $elem       = $1;
      $massnumber = $2;

      die "ERROR: $whoami: no data found for element '$elem'"    unless defined $self->{$elem};
      die "ERROR: $whoami: no data found for isotope '$isotope'" unless defined $self->{$elem}{isotopes}{$isotope};
   }
   else
   {
      die "ERROR: $whoami: invalid syntax of isotope '$isotope'";
   }

   return $massnumber;
}
################################################################################


################################################################################
################################################################################
sub get_element_for_isotope
{
   my $whoami  = _whoami();
   my $self    = shift;
   my $isotope = shift;
   my $elem;
   my $massnumber;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);

   if( $isotope =~ /^([A-Z][a-z]*)(\d+)/ )
   {
      $elem = $1;
      $massnumber = $2;

      die "ERROR: $whoami: no data found for element '$elem'" unless defined $self->{$elem};
   }
   else
   {
      die "ERROR: $whoami: invalid syntax of isotope '$isotope'";
   }

   return $elem;
}
################################################################################


################################################################################
################################################################################
sub get_relative_intensity
{
   my $whoami  = _whoami();
   my $self    = shift;
   my $isotope = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);

   my $elem = $self->get_element_for_isotope($isotope);

   return $self->{$elem}{isotopes}{$isotope}{relative_intensities};
}
################################################################################


################################################################################
################################################################################
sub get_number_of_isotopes
{
   my $whoami = _whoami();
   my $self   = shift;
   my $elem   = shift;

   die "ERROR: $whoami: was not called for a " . __PACKAGE__ . " object" unless ref($self) && $self->isa(__PACKAGE__);

   die "ERROR: $whoami: couldn't find isotope data for element '$elem', check your data source" unless defined $self->{$elem};

   return $self->{$elem}{num_isotopes};
}
################################################################################


################################################################################
################################################################################
sub _process_data_sets
{
   my $whoami = _whoami();
   my $ds    = shift;
   my $self  = shift;

   $ds =~ s/^\s*//;
   $ds =~ s/\s*$//;
   my ($str_iso_name, $str_nat_abu) = split ':', $ds;

   # trim strings
   $str_iso_name =~ s/^\s*//;
   $str_iso_name =~ s/\s*$//;
   $str_nat_abu  =~ s/^\s*//;
   $str_nat_abu  =~ s/\s*$//;
   #foreach my $field (@data_fields)
   #{
   #   print "field: $field\n";
   #}

   die "ERROR: $whoami: syntax error in isotope data '$ds'" unless $str_iso_name && $str_nat_abu;

   ############################################################################
   # we expect 7 fields
   # 1. name of isotopes
   ############################################################################
   my @iso_names = split /\s+/, $str_iso_name;
   my $num_isotopes = @iso_names;

   die "ERROR: $whoami: name of isotope not defined" unless defined $iso_names[0];
   ############################################################################

   ############################################################################
   # get element name
   ############################################################################
   my $elem_name;
   my $iso_massnumber;
   if( $iso_names[0] =~ /^([A-Z][a-z]*)(\d+)$/ )
   {
      $elem_name = $1;
      $iso_massnumber  = $2;
      # warn "INFO: $whoami: found isotopes $1$2 for element '$elem_name'\n";
   }
   else
   {
      die "ERROR: $whoami: isotope name '$iso_names[0]' is invalid";
   }

   $self->{$elem_name}{isotopes}{$iso_names[0]}{massnumber} = $iso_massnumber;
   ############################################################################

   ############################################################################
   # extract info from isotope names
   ############################################################################
   for( my $i = 1; $i < $num_isotopes; $i++ )
   {
      if( $iso_names[$i] =~ /^([A-Z][a-z]*)(\d+)$/ )
      {
         my $e_name = $1;
         my $i_massnumber = $2;
         if( $elem_name ne $1 )
         {
            die "ERROR: $whoami: found isotope for element '$e_name' but expected '$elem_name'";
         }
         $self->{$elem_name}{isotopes}{$iso_names[$i]}{massnumber} = $i_massnumber;
      }
      else
      {
         die "ERROR: $whoami: isotope name '$iso_names[$i]' is invalid";
      }
   }
   $self->{$elem_name}{num_isotopes} = $num_isotopes;
   ############################################################################
   

   ############################################################################
   # 7. relative intensity (most abundant is 100%)
   ############################################################################
   my @rel_intensities = split /\s+/, $str_nat_abu;
   my $num_rel_intensities = @rel_intensities;
   if( $num_isotopes != $num_rel_intensities )
   {
      warn  "ERROR: $whoami: isotope names: @iso_names";
      warn  "ERROR: $whoami: number of protons: @rel_intensities";
      die "ERROR: $whoami: ($iso_names[0]) differing number of isotopes ($num_isotopes) and number of abundance fields ($num_rel_intensities)";
   }

   #if( $rel_intensities[0] != 100 )
   #{
   #   warn "ERROR: $whoami: relative intensity of first entry is invalid ($rel_intensities[0]),\n";
   #   warn "                this value must always by 100.0 and the others relative to it,\n";
   #   warn "                even it this would mean values larger than 100.0\n";
   #   die  "                execution aborted.\n";
   #}

   # copy data to hash/object variable
   $self->{$elem_name}{lowest_isotope} = $iso_names[0];

   my $sum = 0;
   for( my $i = 0; $i < $num_rel_intensities; $i++ )
   {
      $self->{$elem_name}{isotopes}{$iso_names[$i]}{relative_intensities} = $rel_intensities[$i];
      if( $rel_intensities[$i] < 0.0 )
      {
         die "ERROR: $whoami: relative intensity ($rel_intensities[$i]) for isotope '$iso_names[$i]' is invalid.\n";
      }
      $sum += $rel_intensities[$i];
   }
   ############################################################################

   ############################################################################
   die "ERROR: $whoami: sum of natural abundance for element ($elem_name) is not equal to 1 ($sum)\n" if abs($sum - 1) > 1e-8;
   ############################################################################

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

__DATA__
C12 C13: 0.989299998329072 0.0107000016709277
H1 H2: 0.999885003225779 0.000114996774220997
N14 N15: 0.996359998019236 0.00364000198076368
O16 O17 O18: 0.997569999287336 0.000380004339828525 0.00204999637283548
Si28 Si29 Si30: 0.922229995288972 0.0468500031000761 0.0309200016109514
S32 S33 S34 S36: 0.949900003551676 0.00750000197104251 0.0424999985039075 9.99959733738849e-05
P31: 1.0
