## Process this file with automake to produce Makefile.in

fieldsdir   = $(pkgdatadir)/fields
fields_DATA = # Data installed by 'make check'
noinst_DATA = # Data kept with distribution but not installed

# Steady state laminar channel flow fields at a variety of orders
fields_DATA += channel_k06.h5
fields_DATA += channel_k07.h5
fields_DATA += channel_k08.h5
fields_DATA += channel_k09.h5

# The zero-zero modes from Coleman JCP 1995-like turbulent channels.
# Useful for testing implicit solver or relaxation time to laminar profile.
fields_DATA += coleman3k01.00.h5 coleman3k05.00.h5
fields_DATA += coleman3k15.00.h5 coleman3k30.00.h5
fields_DATA += coleman5k01.00.h5 coleman5k05.00.h5
fields_DATA += coleman5k15.00.h5 coleman5k30.00.h5

# Temporally homogenized boundary layers at Coleman-like conditions
fields_DATA += temporal_1e-1_k08.h5
fields_DATA += temporal_1e-2_k08.h5
fields_DATA += temporal_1e-3_k08.h5

# Spatiotemporally homogenized boundary layers at Coleman-like conditions
fields_DATA += spatiotemporal_1e-1_k06.h5                      # ZPG
fields_DATA += spatiotemporal_consistent_1e-1_k06.h5           # ZPG
fields_DATA += spatiotemporal_consistent_1e-1_pexi-1e-2_k06.h5 # FPG
fields_DATA += spatiotemporal_consistent_5e-2_pexi-1e-2_k06.h5 # FPG

# Coarse-but-turbulent fields suitable for development or production spin-up
# See Redmine #3058 if you need the backstory on obtaining these fields.
# Beware these employ bulk viscosity ('/alpha') to permit coarse grids!
fields_DATA += coarse3D.channel.h5
fields_DATA += coarse3D.couette.h5
fields_DATA += coarse3D.couette.forced.h5
fields_DATA += coarse3D.bl.isothermal.zpg.h5

# Coarse-but-turbulent fields per CEV ISS laminar edge (4.134m) conditions
# See Redmine #3088 if you need the backstory on obtaining these fields.
# Beware these employ bulk viscosity ('/alpha') to permit coarse grids!
# fields_DATA += coarse3D.cevisslamedge.isothermal.zpg.h5    # TODO
# fields_DATA += coarse3D.cevisslamedge.isothermal.fpg.h5    # TODO
# fields_DATA += coarse3D.cevisslamedge.cold.fpg.h5          # TODO
# fields_DATA += coarse3D.cevisslamedge.cold.blowing.fpg.h5  # TODO

# Mean profiles for some of coarse3D.* so one might see what one is using.
# Using postproc/gplot is recommended to interrogate this data.  For example,
#     gplot -c using '"y":"bar_u"' with linespoints ::: coarse3D.*.mean.gz
noinst_DATA += coarse3D.channel.mean.gz
noinst_DATA += coarse3D.couette.mean.gz
noinst_DATA += coarse3D.couette.forced.mean.gz
noinst_DATA += coarse3D.bl.isothermal.zpg.mean.gz

# Legacy fields for checking backwards-compatibility of restarting
noinst_DATA += legacy_r22804.h5
noinst_DATA += legacy_r38808.h5

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
