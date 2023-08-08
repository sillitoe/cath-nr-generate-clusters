
use 5.006;
use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME                => 'cath-nr-generate-clusters',
    AUTHOR              => q{Ian Sillitoe <i.sillitoe@ucl.ac.uk>},
    VERSION             => '0.0.1',
    ABSTRACT            => 'Generate non-redundant sequence clusters from BLAST database',
    ($ExtUtils::MakeMaker::VERSION >= 6.3002
      ? ('LICENSE'=> 'perl')
      : ()),
    PREREQ_PM => {
      'Data::Dumper',
      'English',
      'feature',
      'List::MoreUtils',
      'Path::Class',
      'Getopt::Long',
    },
);
