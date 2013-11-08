## Process this file with automake to produce Makefile.in

fieldsdir   = $(pkgdatadir)/fields
fields_DATA = # Data installed by 'make check'
noinst_DATA = # Data kept with distribution but not installed

# Steady state laminar channel flow fields at a variety of orders
fields_DATA += channel_k06.h5 channel_k07.h5 channel_k08.h5 channel_k09.h5

# The zero-zero modes from Coleman JCP 1995-like turbulent channels.
# Useful for testing implicit solver or relaxation time to laminar profile.
noinst_DATA += coleman3k01.00.h5 coleman3k05.00.h5
noinst_DATA += coleman3k15.00.h5 coleman3k30.00.h5
noinst_DATA += coleman5k01.00.h5 coleman5k05.00.h5
noinst_DATA += coleman5k15.00.h5 coleman5k30.00.h5

# Temporally homogenized boundary layers at Coleman-like conditions
noinst_DATA += temporal_1e-1_k08.h5 temporal_1e-2_k08.h5 temporal_1e-3_k08.h5

# Spatiotemporally homogenized boundary layers at Coleman-like conditions
noinst_DATA += spatiotemporal_1e-1_k06.h5 spatiotemporal_consistent_1e-1_k06.h5

# Legacy fields for checking backwards-compatibility of restarting
noinst_DATA += legacy_r22804.h5
noinst_DATA += legacy_r38808.h5

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
