#!/usr/bin/perl

use strict;
use warnings;

use Test::More;
use FindBin;
use Data::Tabix;

my $fn = 'test.bg.gz';

chdir $FindBin::Bin;

require_ok ("Data::Tabix");

ok( my $r = Data::Tabix->new($fn), "loaded input" );
ok( my @names = $r->names, "names()" );
ok( scalar @names == 2, "name count" );

ok( my $res = $r->query('seq2', 112, 154), "query()" );
ok( scalar @$res  == 1, "query count" );
ok( $res->[0]->[3] == 0, "depth check" );

ok( $res = $r->query('NC_019471.2', 446617, 446632), "query()" );
ok( scalar @$res  == 2, "query count" );
ok( $res->[1]->[3] == 7, "depth check" );

ok( ! defined $r->query('foobar', 1, 4), "bad query()" );

ok( $res = $r->query('seq2'), "query w/ no coords" );
ok( scalar @$res  == 6, "query count" );

done_testing();
