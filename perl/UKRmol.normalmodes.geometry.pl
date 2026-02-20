# Geometries can be specified either manually in the array 'geometries'
# by copying what is between start copy and end copy for each geometry
# or automatically as shown below the hash array %geometry
#
# Note that 'length_unit' should be set to the actual length unit used for values in 'atoms' array,
# later (in &generate_geometries) $model{'r_unit'} is set also to this unit for consistency reasons
#
# Note that there's a restriction on the orientation of the molecules in the scripts (although there isn't one in the codes themselves)
# Your molecular target should be oriented as follows for the corresponding point groups ('any orientation' still requires the axes to
# be along the cartesian axes):
# D2h: any orientation
# C2v: C2 axis along Z
# C2h: C2 axis along Z
# D2:  any orientation
# C2:  C2 axis along Z
# Cs:  Molecule in the YZ plane
# Ci:  any orientation
# C1:  any orientation

# ----------- CHANGE THESE, make sure that the geometries agree with the above restrictions
my $SUFF = "ASYMSTRETCH";
my $CORR_CM = 0; # 0: no com correction. 1: com correction
# my $INFILE = "/home/josh/Workspace/ResearchCodes/UKRMol+/Input/H2OX-e/vib/bend/H2OX.bend.xyz";
my $INFILE = "H2OX.asymstretch.xyz";
# ------------------------------

# -- Read the MOL.MODE.XYZ file to get the displacements for this normal mode
use strict;
use warnings;

open(my $fh, "<", $INFILE) or die "Can't open $INFILE: $!";

local $/ = "\n"; # <-- local record separator override IMPORTANT

my ($natoms, $ndisp, $units);
my @geoms;

# -- trim space
sub trim {
  my ($s) = @_;
  $s =~ s/^\s+|\s+$//g;
  return $s;
}

# -- read nonempty line
sub next_nonempty_line {
  my ($fh) = @_;
  while (my $l = <$fh>) {
    $l = trim($l);
    next if $l eq "";
    return $l;
  }
  return undef;
}

while (my $line = next_nonempty_line($fh)) {

  # -- header fields (can appear anywhere before blocks)
  $natoms = 0 + $1 if $line =~ /^NATOMS\s*=\s*(\d+)/i;
  $ndisp  = 0 + $1 if $line =~ /^NDISP\s*=\s*(\d+)/i;
  $units  = lc($1) if $line =~ /^XYZ\s+units\s*=\s*([A-Za-z0-9_]+)/i;

  # -- Q block
  next unless $line =~ /^Q\s*=\s*([+-]?\d+(?:\.\d+)?)/i;
  die "NATOMS not set before first Q block" unless defined $natoms;

  my $Q = 0.0 + $1;
  my @atoms;

  while (@atoms < $natoms) {
    my $a = next_nonempty_line($fh);
    die "EOF while reading atoms for Q=$Q" unless defined $a;

    my ($sym, $x, $y, $z) = split(/\s+/, $a);
    die "Bad atom line: '$a'" unless defined $z;

    push @atoms, [ $sym, 0.0+$x, 0.0+$y, 0.0+$z ];
  }

  push @geoms, {
    description  => sprintf("%+7.3f", $Q),
    gnuplot_desc => "Q=".sprintf("%+.3f", $Q),
    atoms        => \@atoms,
  };
}

close($fh);

# units: 1=angstrom else bohr(0)
my $length_unit = (defined($units) && $units =~ /^ang/i) ? 1 : 0;

our %geometry;
%geometry = (
  suffix            => $SUFF,
  correct_cm        => $CORR_CM,
  length_unit       => $length_unit,
  geometry_labels   => "    Q     ",
  geometries        => \@geoms,
  start_at_geometry => 1,
  stop_at_geometry  => 0,
);

warn sprintf("Parsed %d geometries (NDISP=%s) from %s\n",
  scalar(@geoms), (defined $ndisp ? $ndisp : "NA"), $INFILE);
