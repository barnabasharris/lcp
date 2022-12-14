module purge
module load beta-modules
module load r/r-4.1.1_bc-3.13
module load nano/2.4.2
module load gerun

export TZ="Europe/London"
export TMPDIR="/home/tcrnbgh/Scratch/tmp"

echo "TMPDIR=/home/tcrnbgh/Scratch/tmp" > /home/tcrnbgh/Scratch/lcp/.Renviron
cp /home/tcrnbgh/Scratch/lcp/.Renviron /home/tcrnbgh/.Renviron




