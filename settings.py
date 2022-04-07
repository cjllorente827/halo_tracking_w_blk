HOME = "/mnt/home/llorente"
SCRATCH = "/mnt/gs18/scratch/users/llorente"
ENZODIR = f"{SCRATCH}/bigbox/25Mpc_256"
SCRATCHCACHE = f"{SCRATCH}/halo_tracking"
HOMECACHE = f"{HOME}/halo_tracking_pipeline/cache"
CATALOGDIR = f"{HOME}/25Mpc_256_halo_prospecting"
REDSHIFTS = [2,3,4,5,6]

# Indexed by redshift
ENZO_OUTPUTS = {
    2: f"{ENZODIR}/RD0111/RD0111",
    3: f"{ENZODIR}/RD0077/RD0077",
    4: f"{ENZODIR}/RD0053/RD0053",
    5: f"{ENZODIR}/RD0036/RD0036",
    6: f"{ENZODIR}/RD0022/RD0022",
}

# Indexed by redshift
HALO_CATALOGS = {
    2 : f"{CATALOGDIR}/rockstar_halos/out_4.list",
    3 : f"{CATALOGDIR}/rockstar_halos/out_3.list",
    4 : f"{CATALOGDIR}/rockstar_halos/out_2.list",
    5 : f"{CATALOGDIR}/rockstar_halos/out_1.list",
    6 : f"{CATALOGDIR}/rockstar_halos/out_0.list"
}

ANNAS_CATALOG = f"{CATALOGDIR}/halo_catalog_sorted.ascii"