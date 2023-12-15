#!/bin/bash
#-------- Set up and run real for a WRF met-only run --------
#
# Louis Marelle, 2023/12/15
#

# Resources used
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00


#-------- Input --------
CASENAME='WRF_MET_SIMPLE'
CASENAME_COMMENT=''

# Root directory with the compiled WRF executables (main/wrf.exe and main/real.exe)
WRFDIR=~/WRF/src/WRF-Chem-Polar/WRFV4
WRFVERSION='met'

# Simulation start year and month
yys=2012
mms=01
dds=01
hhs=00
# Simulation end year, month, day, hour
yye=2012
mme=01
dde=05
hhe=00

NAMELIST="namelist.input.YYYY"


#-------- Parameters --------
# Root directory for WRF input/output
OUTDIR_ROOT="/data/$(whoami)/WRF/WRF_OUTPUT"
SCRATCH_ROOT="/scratchu/$(whoami)"


#-------- Set up job environment --------
# Load modules used for WRF compilation
module purge
module load gcc/11.2.0
module load openmpi
module load netcdf-c/4.7.4
module load netcdf-fortran/4.5.3
module load hdf5/1.10.7
module load jasper

# Set run start and end date
date_s="$yys-$mms-$dds"
date_e="$yye-$mme-$dde"


#-------- Set up real input and output directories & files  --------
# Run id
ID="$(date +"%Y%m%d").$SLURM_JOBID"

# Case name for the output folder
if [ -n "$CASENAME_COMMENT" ]; then 
  CASENAME_COMMENT="_${CASENAME_COMMENT}"
fi

# Directory containing real.exe output (e.g. wrfinput_d01, wrfbdy_d01 files)
REALDIR="${OUTDIR_ROOT}/real_${CASENAME}${CASENAME_COMMENT}_$(date -d "$date_s" "+%Y")"
mkdir "$REALDIR"
WPSDIR="${OUTDIR_ROOT}/met_em_${CASENAME}_$(date -d "$date_s" "+%Y")"

# Also create a temporary scratch run directory
SCRATCH="$SCRATCH_ROOT/real_${CASENAME}${CASENAME_COMMENT}_$(date -d "$date_s" "+%Y").${ID}.scratch"
rm -rf "$SCRATCH"
mkdir $SCRATCH
cd $SCRATCH

# Init spectral nudging parameters - we only nudge
# the 1000 km = 1000000m scale
nudging_scale=1000000
wrf_dx=$(sed -n -e 's/^[ ]*dx[ ]*=[ ]*//p' "$SLURM_SUBMIT_DIR/${NAMELIST}" | sed -n -e 's/,.*//p')
wrf_dy=$(sed -n -e 's/^[ ]*dy[ ]*=[ ]*//p' "$SLURM_SUBMIT_DIR/${NAMELIST}" | sed -n -e 's/,.*//p')
wrf_e_we=$(sed -n -e 's/^[ ]*e_we[ ]*=[ ]*//p' "$SLURM_SUBMIT_DIR/${NAMELIST}" | sed -n -e 's/,.*//p')
wrf_e_sn=$(sed -n -e 's/^[ ]*e_sn[ ]*=[ ]*//p' "$SLURM_SUBMIT_DIR/${NAMELIST}" | sed -n -e 's/,.*//p')
xwavenum=$(( (wrf_dx * wrf_e_we) / $nudging_scale))
ywavenum=$(( (wrf_dy * wrf_e_sn) / $nudging_scale))

# Write the info on input/output directories to run log file
echo "Running real.exe from $WRFDIR"
echo "Running on scratchdir $SCRATCH"
echo "Writing output to $REALDIR"
echo "WPS dir $WPSDIR"
echo "Running from $date_s to $date_e"

# Save this slurm script to the output directory
cp $0 "$REALDIR/jobscript_real.sh"


#-------- Run real --------
cd $SCRATCH

#---- Copy all needed files to scrach space
# Input files from run setup directory
cp "$SLURM_SUBMIT_DIR/"* "$SCRATCH/"
# Executables and WRF aux files from WRFDIR
cp "$WRFDIR/run/"* "$SCRATCH/"
cp "$WRFDIR/../executables/real.exe.$WRFVERSION" "$SCRATCH/real.exe"
# met_em WPS files from WPSDIR
cp "${WPSDIR}/met_em.d"* "$SCRATCH/"


#---- Run real.exe
# Prepare the real.exe namelist, set up run start and end dates
cp "$SLURM_SUBMIT_DIR/${NAMELIST}" namelist.input
sed -i "s/__STARTYEAR__/${yys}/g" namelist.input
sed -i "s/__STARTMONTH__/${mms}/g" namelist.input
sed -i "s/__STARTDAY__/${dds}/g" namelist.input
sed -i "s/__STARTHOUR__/${hhs}/g" namelist.input
sed -i "s/__ENDYEAR__/${yye}/g" namelist.input
sed -i "s/__ENDMONTH__/${mme}/g" namelist.input
sed -i "s/__ENDDAY__/${dde}/g" namelist.input
sed -i "s/__ENDHOUR__/${hhe}/g" namelist.input
sed -i "s/__XWAVENUM__/$xwavenum/g" namelist.input
sed -i "s/__YWAVENUM__/$ywavenum/g" namelist.input
echo " "
echo "-------- jobscript: run real.exe --------"
echo " "
mpirun ./real.exe
# Check the end of the log file in case real crashes
tail -n20 rsl.error.0000


#-------- Transfer data  --------
# Clean up
rm -f met_em*
# Transfer files to the output dir
cp rsl* "$REALDIR/"
cp *d0* "$REALDIR/"
cp namelist.input "$REALDIR/namelist.input.real"

# Remove scratch dir
rm -rf "$SCRATCH"

