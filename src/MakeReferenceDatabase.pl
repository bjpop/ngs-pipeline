#!/usr/bin/env perl

use strict;
use warnings;
use IndexReference;
use Options;
use Timer;

# Process command line arguments. Remove any arguments from command line vector.
my (%OPTIONS, $REFERENCE, @SEQUENCES) = Options::get_options();

# Step 1. Create reference database.
IndexReference::mkRefDataBase($REFERENCE, "-a $OPTIONS{i}");

Timer::timer('Reference database created');

# Step 2. Index the reference.
IndexReference::indexReference($REFERENCE);

Timer::timer('Reference database indexed');
