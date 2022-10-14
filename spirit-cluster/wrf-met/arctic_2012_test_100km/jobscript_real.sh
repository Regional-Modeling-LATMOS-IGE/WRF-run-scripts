#!/bin/bash
#-------- Set up and run real for a WRF met-only run --------
#
# Louis Marelle, 2022/10/04
#

# Resources used
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00


#-------- Input --------
CASENAME='WRF_MET_SIMPLE'

# Root directory with the compiled WRF executables (main/wrf.exe and main/real.exe)
WRFDIR=~/WRF/src/WRF-Chem-Polar/WRFV4

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
# Root directory for WRF output
OUTDIR_ROOT="/data/$(whoami)/WRF/WRF_OUTPUT"
SCRATCH_ROOT="/scratchu/$(whoami)"


#-------- Set up job environment --------
# Load modules used for WRF compilation
module purge
module load gcc/11.2.0
module load openmpi/4.0.7 
module load netcdf-c/4.7.4
module load netcdf-fortran/4.5.3
module load hdf5/1.10.7
module load jasper

# Set run start and end date
date_s="$yys-$mms-$dds"
date_e="$yye-$mme-$dde"


# -------- Sanity checks on inputs --------
echo ""
# All of these inputs should be integers
if ! [ "$yys" -eq "$yys" ] || ! [ "$mms" -eq "$mms" ] ||  ! [ "$dds" -eq "$dds" ] \
|| ! [ "$yye" -eq "$yye" ] || ! [ "$mme" -eq "$mme" ] ||  ! [ "$dde" -eq "$dde" ]; then
  echo "Error, inputs to this script should be integers" >&2
  exit 1
fi
# Dates should be in the YYYY-MM-DD format
if (( ${#yys} !=4 | ${#mms} !=2 | ${#dds} !=2  )); then
  echo "Error, start year, month, date format must be YYYY, MM, DD; now $yys, $mms, $dds" >&2
  exit 1
fi
if (( ${#yye} !=4 | ${#mme} !=2 | ${#dde} !=2  )); then
  echo "Error, end year, month, date format must be YYYY, MM, DD; now $yye, $mme, $dde" >&2
  exit 1
fi
# Start date should be before end date
if  (( $(date -d "$date_s" "+%s") >= $(date -d "$date_e" "+%s") )); then
  echo "Error: start date $date_s >= end date $date_e" >&2
  exit 1
fi


#-------- Set up real input and output directories & files  --------
# Run id
ID="$(date +"%Y%m%d").$SLURM_JOBID"

# Directory containing real output (e.g. wrfinput_d01, wrfbdy_d01 files)
REALDIR="${OUTDIR_ROOT}/real_${CASENAME}_$(date -d "$date_s" "+%Y")"
mkdir "$REALDIR"
WPSDIR="${OUTDIR_ROOT}/met_em_${CASENAME}_$(date -d "$date_s" "+%Y")"

# Also create a temporary run directory
SCRATCH="$SCRATCH_ROOT/real_${CASENAME}_$(date -d "$date_s" "+%Y").scratch"
rm -rf "$SCRATCH"
mkdir $SCRATCH
cd $SCRATCH

# Write the info on input/output directories to run log file
echo "Running real.exe from $WRFDIR"
echo "Running on scratchdir $SCRATCH"
echo "Writing output to $REALDIR"
echo "WPS dir $WPSDIR"
echo "Running from $date_s to $date_e"

# Save this slurm script to the output directory
cp $0 "$REALDIR/jobscript_real.sh"


#-------- Run real --------
# Copy the WRF run directory (contains auxilliary files etc.) and the real.exe
# executable to $SCRATCH/
cp "$SLURM_SUBMIT_DIR/"* "$SCRATCH/"
# Executables and WRF aux files from WRFDIR
cp "$WRFDIR/run/"* "$SCRATCH/"
cp "$WRFDIR/main/real.exe" "$SCRATCH/"
# We run the version real.exe.$WRFVERSION in $WRFDIR/../executables
# WRFVERSION='met'
# cp "$WRFDIR/../executables/real.exe.$WRFVERSION" "$SCRATCH/real.exe"

#  Copy and prepare the WRF namelist, set up run start and end dates
cp "$SLURM_SUBMIT_DIR/${NAMELIST}" namelist.input
sed -i "s/__STARTYEAR__/${yys}/g" namelist.input
sed -i "s/__STARTMONTH__/${mms}/g" namelist.input
sed -i "s/__STARTDAY__/${dds}/g" namelist.input
sed -i "s/__STARTHOUR__/${hhs}/g" namelist.input
sed -i "s/__ENDYEAR__/${yye}/g" namelist.input
sed -i "s/__ENDMONTH__/${mme}/g" namelist.input
sed -i "s/__ENDDAY__/${dde}/g" namelist.input
sed -i "s/__ENDHOUR__/${hhe}/g" namelist.input

# Copy the input files from WPS
date_s_met=$(date +"%Y%m%d" -d "$date_s")
while (( $(date -d "$date_s_met" "+%s") <= $(date -d "$date_e" "+%s") )); do
  yys_met=${date_s_met:0:4}
  mms_met=${date_s_met:4:2}
  dds_met=${date_s_met:6:2}
  cp "${WPSDIR}/met_em.d0"*"$yys_met-$mms_met-$dds_met"* "$SCRATCH/"
  date_s_met=$(date +"%Y%m%d" -d "$date_s_met + 1 day")
done

# Run real.exe --------
echo " "
echo "-------- jobscript: run real.exe --------"
echo " "
mpirun ./real.exe

# Check the end of the log file in case the code crashes
tail -n20 rsl.error.0000

# Clean up
rm -f met_em*

# Transfer files to the output dir
cp rsl* "$REALDIR/"
cp *d0* "$REALDIR/"
cp namelist.input "$REALDIR/namelist.input.real"

# Clean up
rm -rf "$SCRATCH"

