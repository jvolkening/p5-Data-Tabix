use 5.012;

use strict;
use warnings FATAL => 'all';

use Config;
use ExtUtils::MakeMaker;

# require 64-bit integer support
if (! $Config{use64bitint}) {
    warn "64-bit integer support required for installation\n";
    exit(0);
}

WriteMakefile(
    NAME             => 'Data::Tabix',
    AUTHOR           => q{Jeremy Volkening <jdv@base2bio.com>},
    VERSION_FROM     => 'lib/Data/Tabix.pm',
    ABSTRACT_FROM    => 'lib/Data/Tabix.pm',
    LICENSE          => 'GPL_3',
    PL_FILES         => {},
    MIN_PERL_VERSION => 5.012,
    CONFIGURE_REQUIRES => {
        'Config'              => 0,
        'ExtUtils::MakeMaker' => 0,
    },
    BUILD_REQUIRES => {
        'Test::More' => 0,
    },
    PREREQ_PM => {
        'Compress::BGZF' => 0.004,
    },
    EXE_FILES => [qw|
    |],
    META_MERGE => {
        'meta-spec' => { version => 2 },
        resources => {
            repository => {
                type => 'git',
                url => 'https://github.com/jvolkening/p5-Data-Tabix.git',
                web => 'https://github.com/jvolkening/p5-Data-Tabix',
            },
        },
    },
    dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean => { FILES => 'Data-Tabix-*' },
);
