#!/bin/bash
# =====================================================================
# Basado en BITACORA v1.4 - runBITACORA.sh (Joel Vizueta, github/jvizueta@ub.edu).
# Unica modificacion: la seccion de configuracion se lee de VARIABLES DE ENTORNO
# (que Snakemake exporta por muestra). La seccion "PIPELINE - CODE" es identica
# al upstream. Si actualizas BITACORA, re-sincroniza esa seccion.
# =====================================================================
VERSION=1.4

##########################################################
##       CONFIGURACION DESDE VARIABLES DE ENTORNO       ##
##########################################################
# Perl debe estar instalado; BLAST+ y HMMER vienen del env conda (--use-conda).

SCRIPTDIR="${SCRIPTDIR:?Falta SCRIPTDIR}"        # carpeta Scripts de BITACORA
GEMOMAP="${GEMOMAP:-}"                            # jar de GeMoMa (solo si GEMOMA=T)

NAME="${NAME:?Falta NAME}"
GENOME="${GENOME:?Falta GENOME}"                 # ruta absoluta
GFFFILE="${GFFFILE:?Falta GFFFILE}"              # ruta absoluta
PROTFILE="${PROTFILE:?Falta PROTFILE}"           # ruta absoluta
QUERYDIR="${QUERYDIR:?Falta QUERYDIR}"           # carpeta con las FPDB (por familia)

CLEAN="${CLEAN:-T}"
EVALUE="${EVALUE:-1e-5}"
THREADS="${THREADS:-4}"
GEMOMA="${GEMOMA:-F}"
MAXINTRON="${MAXINTRON:-0}"
GENOMICBLASTP="${GENOMICBLASTP:-F}"
ADDFILTER="${ADDFILTER:-T}"
FILTERLENGTH="${FILTERLENGTH:-60}"
RETAINNONFILTER="${RETAINNONFILTER:-F}"

##########################################################
##            PIPELINE - CODE (identico v1.4)           ##
##########################################################

echo -e "\n#######################  Running BITACORA  #######################";
echo "BITACORA version $VERSION";
date

# Checking if provided data is ok

if [[ ! -f $SCRIPTDIR/check_data.pl ]] ; then
	echo -e "BITACORA can't find Scripts folder in $SCRIPTDIR. Be sure to add also Scripts at the end of the path as /path/Scripts";
	echo -e "BITACORA died with error\n";
	exit 1;
fi

if [ $GEMOMA == "T" ] ; then
	perl $SCRIPTDIR/check_data.pl $GFFFILE $GENOME $PROTFILE $QUERYDIR $GEMOMA $GEMOMAP 2>BITACORAstd.err
fi

if [ $GEMOMA != "T" ] ; then
	perl $SCRIPTDIR/check_data.pl $GFFFILE $GENOME $PROTFILE $QUERYDIR $GEMOMA 2>BITACORAstd.err
fi

ERRORCHECK="$(grep -c 'ERROR' BITACORAstd.err)"

if [ $ERRORCHECK != 0 ]; then
	cat BITACORAstd.err;
	echo -e "BITACORA died with error\n";
	exit 1;
fi


# Run step 1

perl $SCRIPTDIR/runanalysis.pl $NAME $PROTFILE $QUERYDIR $GFFFILE $GENOME $EVALUE $THREADS 2>>BITACORAstd.err

ERRORCHECK="$(grep -c 'ERROR' BITACORAstd.err)"

if [ $ERRORCHECK != 0 ]; then
	cat BITACORAstd.err;
	echo -e "BITACORA died with error\n";
	exit 1;
fi


# Run step 2

if [ $GEMOMA == "T" ] ; then
	if [ $GENOMICBLASTP == "T" ] ; then
		perl $SCRIPTDIR/runanalysis_2ndround_v2_genomic_withgff_gemoma.pl $NAME $PROTFILE $QUERYDIR $GENOME $GFFFILE $EVALUE $MAXINTRON $THREADS $GEMOMAP 2>>BITACORAstd.err 2>BITACORAstd.err
	fi

	if [ $GENOMICBLASTP != "T" ] ; then
		perl $SCRIPTDIR/runanalysis_2ndround_genomic_withgff_gemoma.pl $NAME $PROTFILE $QUERYDIR $GENOME $GFFFILE $EVALUE $MAXINTRON $THREADS $GEMOMAP 2>>BITACORAstd.err 2>BITACORAstd.err
	fi
fi

if [ $GEMOMA != "T" ] ; then
	if [ $GENOMICBLASTP == "T" ] ; then
		perl $SCRIPTDIR/runanalysis_2ndround_v2_genomic_withgff.pl $NAME $PROTFILE $QUERYDIR $GENOME $GFFFILE $EVALUE $MAXINTRON $THREADS 2>>BITACORAstd.err 2>BITACORAstd.err
	fi

	if [ $GENOMICBLASTP != "T" ] ; then
		perl $SCRIPTDIR/runanalysis_2ndround_genomic_withgff.pl $NAME $PROTFILE $QUERYDIR $GENOME $GFFFILE $EVALUE $MAXINTRON $THREADS 2>>BITACORAstd.err 2>BITACORAstd.err
	fi
fi

ERRORCHECK="$(grep -c 'ERROR' BITACORAstd.err)"

if [ $ERRORCHECK != 0 ]; then
	cat BITACORAstd.err;
	echo -e "BITACORA died with error\n";
	exit 1;
fi

ERRORCHECK="$(grep -c 'Segmentation' BITACORAstd.err)"

if [ $ERRORCHECK != 0 ]; then
	cat BITACORAstd.err;
	echo -e "BITACORA died with error\n";
	exit 1;
fi


# Run additional filtering and clustering

if [ $ADDFILTER == "T" ] ; then
	perl $SCRIPTDIR/runfiltering.pl $NAME $QUERYDIR $FILTERLENGTH 2>>BITACORAstd.err 2>BITACORAstd.err
fi

ERRORCHECK="$(grep -c 'ERROR' BITACORAstd.err)"

if [ $ERRORCHECK != 0 ]; then
	cat BITACORAstd.err;
	echo -e "BITACORA died with error\n";
	exit 1;
fi


# Cleaning

if [ $RETAINNONFILTER == "T" ]; then
	perl $SCRIPTDIR/runcleaning_allcopies.pl $NAME $QUERYDIR
	echo -e "Cleaning output folders\n";
fi

if [ $RETAINNONFILTER != "T" ]; then
	if [ $CLEAN == "T" ]; then
	perl $SCRIPTDIR/runcleaning.pl $NAME $QUERYDIR
	echo -e "Cleaning output folders\n";
	fi
fi


rm BITACORAstd.err

echo -e "BITACORA completed without errors :)";
date
