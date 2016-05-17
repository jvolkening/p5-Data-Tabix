package Data::Tabix 0.001;

use 5.012;

use strict;
use warnings;

use Carp;
use Compress::BGZF::Reader;
use List::Util qw/any max/;

use constant TABIX_MAGIC => 'TBI' . chr(1);
use constant HEAD_BYTES  => 36; # not including extra fields

sub new {

    #-------------------------------------------------------------------------
    # ARG 0 : input filename (must be BGZF compressed and tabix-indexed)
    #-------------------------------------------------------------------------
    # RET 0 : Data::Tabix object
    #-------------------------------------------------------------------------
    
    my ($class, $fn_in) = @_;
    my $self = bless {}, $class;

    # initialize
    croak "Input file not specified or not found\n"
        if (! -r $fn_in);
    my $fn_idx = $fn_in . '.tbi';
    croak "Index file not found - is the input indexed by tabix?\n"
        if (! -r $fn_idx);

    my $r_in  = Compress::BGZF::Reader->new($fn_in)
        or die "Error opening input for reading: $@\n";
    my $r_idx = Compress::BGZF::Reader->new($fn_idx)
        or die "Error opening index for reading: $@\n";

    $self->{in} = $r_in;
    $self->_load_idx( $r_idx );

    return $self;

}

sub _load_idx {

    my ($self, $r) = @_;

    # load header
    my @h_keys = qw/
        n_ref
        format
        col_seq
        col_beg
        col_end
        meta
        skip
        l_nm
    /;

    my ($magic, @h_vals) = unpack 'A4l<8', $r->read_data( HEAD_BYTES );
    croak "MAGIC mismatch: possibly corrupt Tabix index file\n"
        if ($magic ne TABIX_MAGIC);
    @{$self->{header}}{@h_keys} = @h_vals;

    # load names
    $self->{names} = [ split /\0/, $r->read_data( $self->{header}->{l_nm} ) ];

    # parse indices

    for my $i (0..$self->{header}->{n_ref}-1) {
        my $n_bin = unpack 'l<', $r->read_data(4);
        for my $b (0..$n_bin-1) {
            my ($bin, $n_chunk) = unpack 'L<l<', $r->read_data(8);
            $self->{idx}->{$self->{names}->[$i]}->{$bin}
                = [ unpack 'Q<*', $r->read_data($n_chunk*16) ];
        }
        my $n_intv = unpack 'l<', $r->read_data(4);
        $self->{idx}->{$self->{names}->[$i]}->{lin}
            = [ unpack 'Q<*', $r->read_data($n_intv*8) ];
    }
    $self->{last_pos_field} = max( $self->{header}->{col_beg},
        $self->{header}->{col_end} );

}

sub query {

    my ($self, $name, $start, $end) = @_;

    return undef if (! any {$_ eq $name} @{$self->{names}});
    $start = $start // 1;
    $end   = defined $end && $end < 1<<29 ? $end : 1<<29;
    return undef if ($end < $start);

    my $in = $self->{in}; # BGZF reader

    # use linear index to find lowest expected matching vo
    my $lin_idx = int(($start-1)/2**14);
    my $vo_lowest = $self->{idx}->{$name}->{lin}->[$lin_idx];

    my @matrix;
    my @blocks;

    # find candidate blocks
    BIN:
    for my $bin ( _candidates( $start, $end ) ) {
        next BIN if (! defined $self->{idx}->{$name}->{$bin});
        my @chunks = @{ $self->{idx}->{$name}->{$bin} };
        my $i = 0;
        CHUNK:
        while ($i < $#chunks) {
            my $vo = $chunks[$i++];
            my $f  = $chunks[$i++];
            next CHUNK if ($f < $vo_lowest);
            push @blocks, [$vo, $f];
        }
    }
    return if (! scalar @blocks);

    # merge overlapping blocks
    @blocks = sort {$a->[0] <=> $b->[0]} @blocks;
    my @final = @{ shift @blocks };
    for (@blocks) {
        if ($_->[0] > $final[-1]) {
            push @final, @$_;
        }
        elsif ($_->[1] > $final[-1]) {
            $final[-1] = $_->[1];
        }
    }

    # test for overlap with target region

    # add 1 to start coord for BED-type coordinates
    my $off_start = ($self->{header}->{format} & 0x10000) ? 1 : 0;
    my $i = 0;
    while ($i < $#final) {
        my $vo = $final[$i++];
        my $f  = $final[$i++];
        $in->move_to_vo ($vo);
        while ($vo < $f) {
            my $line = $in->getline;
            chomp $line;
            $vo = $in->get_vo;

            # split was found to be a major slowdown for files with many
            # columns (e.g. some VCFs, etc). Adding the limit here, even with
            # the second split below, prevents problems for these files

            my @fields = split "\t", $line, $self->{last_pos_field}+1;

            # check for overlap
            last if $fields[$self->{header}->{col_beg}-1]+$off_start > $end;
            if ($self->{header}->{col_end} > 0) {
                next if ($fields[$self->{header}->{col_end}-1] < $start);
            }
            else {
                next if ($fields[$self->{header}->{col_beg}-1]+$off_start < $start);
            }
            push @matrix, [split "\t", $line];
        }
    }
                
    return \@matrix;

}

sub names {

    my ($self) = @_;
    return @{ $self->{names} };

}

sub _candidates {

    my ($s, $e) = @_;

    my @bins = (0);
    my $a    =  1;
    my $b    = 26;
    --$e;
    for (0..4) {
        push @bins, ($a+($s>>$b)..$a+($e>>$b));
        $a  = $a*8 + 1;
        $b -= 3;
    }

    return @bins;

}

1;


__END__

=head1 NAME

Data::Tabix - Pure Perl interface to Tabix-indexed flatfiles

=head1 SYNOPSIS

    use Data::Tabix;

    my $reader = Data::Tabix->new( 'some_file.tsv.gz' );
    my $results = $reader->query( $name, $start, $end );
    for my $row (@$results) {
        
        # $row is an arrayref containing the fields of the record
        # (format-dependent)

    }

=head1 DESCRIPTION

C<Data::Tabix> enables querying of Tabix-indexed tab-delimited tables for
region overlap, most commonly used in genomics. Input files must be block gzip
compressed (e.g. using bgzip) and indexed using tabix. Queries specify a
reference name (typically a sequence ID) and optionally 1-based start and end
coordinates. The returned array reference contains a two-dimensional array of
records which overlapped the requested region. For details of Tabix indexing
and retrieval, see the original publication at
L<http://dx.doi.org/10.1093/bioinformatics/btq671>.

This module is written in pure Perl as an alternative to XS-based modules (see below).

=head1 METHODS

=over 4

=item B<new>

    my $reader = Data::Tabix->new( $input_fn );

Creates and returns a new C<Data::Tabix> engine. Takes a file path as the
single required input parameter. The file must be BGZF compressed (e.g. using
bgzip) and indexed using tabix. The index file must be present as
"$input_fn.tbi". During construction the entire index is parsed, which is
moderately slower than the C-based alternatives but is only done once.
Subsequent queries are rapid, such that for a reasonably large number of
queries the different in performance becomes negligible.

=item B<query>

    # returns all records overlapping the region between 1000-2000 bp on chr1
    my $results = $reader->query("chr1", 1000, 2000);
    for (@$results) {
        my $strand = $_->[STRAND_COL];
    }

Search the indexed file and return all records whose extent (based on one or
two coordinate fields) overlaps that requested. Takes a minimum of one
argument (the reference name, typically a sequence ID) and two optional
arguments (starting and ending coordinates of the target region).

Returns a reference to an array with one record per line. Individual records
are also array references containing the fields of the record
(format-specific).

=item B<names>

    my @names = $reader->names;

Returns an array of reference names present in the file. Typically this will
be a list of sequence IDs.

=back

=head1 CAVEATS AND BUGS

This is code is in alpha testing stage and the API is not guaranteed to be
stable.

Please reports bugs to the author.

=head1 SEE ALSO

L<Bio::HTS::Tabix> - XS bindings to the htslib tabix API

=head1 AUTHOR

Jeremy Volkening <jdv@base2bio.com>

=head1 COPYRIGHT AND LICENSE

Copyright 2015-2016 Jeremy Volkening

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

