#!/usr/bin/env bash

usage () {
    cat <<EOF
    Starting from a file of depos, run WCT sim+sigproc+imaging.
    Produce output files and convert them to various vizualization
    formats including paraview and Bee.

    usage:

        test-ssi-viz.sh <detector> <depos.npz> <outdir>/

EOF
    exit 1
}

set -Eeuo pipefail

trap usage SIGINT SIGTERM ERR 

detector="${1:-}"
if [ -z "$detector" ] ; then
    echo "Error: no detector given"
    usage
fi

depos="${2:-}";
if [ -z "$depos" ] ; then
    echo "Error: no input depos file given";
    usage
fi
depos="$(realpath $depos)"
if [ ! -f $depos ] ; then
    echo "Error: no such file: $depos"
    usage
fi

outdir="$(realpath ${3:-.})" 

mydir="$(dirname $(realpath $BASH_SOURCE))"

wctdir="$outdir/$detector/wct"
mkdir -p "$wctdir"

echo "Compile configuration file"
framespat="$wctdir/frames-%(tier)s-%(anode)s.npz"
clustersext="npz"
clusterspat="$wctdir/clusters-%(tier)s-%(anode)s.${clustersext}"
deposout="$wctdir/drifted-depos.npz"

base="$wctdir/depo-ssi-viz"
if [ ! -f "${base}.json" ] ; then
    set -x
    wcsonnet \
        -A detector=$detector \
        -A depofile="$depos" \
        -A frames="$framespat" \
        -A clusters="$clusterspat" \
        $mydir/depo-ssi-viz.jsonnet > "${base}.json" || exit -1
    set +x
fi    
if [ ! -f "${base}.pdf" ] ; then
    set -x
    wirecell-pgraph dotify --no-params "${base}.json" "${base}.pdf" || exit -1
    set +x
    echo "${base}.pdf"
fi

echo "Run wire-cell"
log=$wctdir/log
if [ -f "$log" ] ; then
    echo "wire-cell output found, to rerun, remove:"
    echo "$log"
else
    set -x
    wire-cell \
        -l $log -L debug \
        -c $mydir/depo-ssi-viz.jsonnet \
        -A detector=$detector \
        -A depofile="$depos" \
        -A depos="$deposout" \
        -A frames="$framespat" \
        -A clusters="$clusterspat" || exit -1
    set +x
fi


pltdir="$outdir/$detector/plots"
mkdir -p $pltdir

echo "Plot depos"
for letter in y z
do
    pdf="$pltdir/depos-postdrift-t${letter}.pdf"
    if [ -f "$pdf" ] ; then
        echo "already: $pdf"
    else
        set -x
        wirecell-gen plot-depos -p t${letter}qscat -g0 $deposout $pdf
        set +x
    fi
done

echo "Plot ADC"
for one in $wctdir/frames-adc-*.npz
do
    base="$(basename $one .npz)"
    pdf="$pltdir/${base}.pdf"
    if [ -f "$pdf" ] ; then
        echo "Already have $pdf"
    else
        echo "$one"
        wirecell-util npzls "$one"
        set -x
        wirecell-util npz-to-img \
                      --cmap seismic \
                      --title "$base" \
                      --xtitle 'Relative tick number' \
                      --ytitle 'Relative channel number' \
                      --ztitle 'ADC (median subtracted)' \
                      --vmin '-50' --vmax '50' --mask 0 --dpi 600 \
                      --baseline='median'\
                      -o "$pdf" "$one" || exit -1
        set +x
    fi
done

echo "Plot Signal"
for one in $wctdir/frames-sig-*.npz
do
    base="$(basename $one .npz)"
    pdf="$pltdir/${base}.pdf"
    if [ -f "$pdf" ] ; then
        echo "Already have $pdf"
    else
        echo "$one"
        got="$(wirecell-util npzls "$one")"
        if [ -z "$got" ] ; then
            echo "looks empty: $(ls -l $one)"
            continue
        fi
        set -x
        wirecell-util npz-to-img \
                      --cmap viridis \
                      --title "$base" \
                      --xtitle 'Relative tick number' \
                      --ytitle 'Relative channel number' \
                      --ztitle 'Signal (ionization electrons)' \
                      --vmin '0' --vmax '10000' --mask 0 --dpi 600 \
                      -o "$pdf" "$one" || exit -1
        set +x
    fi
done
         
echo "Make paraview"
pvdir="$outdir/$detector/pv"
mkdir -p $pvdir

echo "Dump drifted blobs to paraview"
deposname="$(basename $depos .npz)"
pvdepos="$pvdir/${deposname}.vtp"
if [ ! -f "$pvdepos" ] ; then
    set -x
    wirecell-img paraview-depos --speed 0 $depos $pvdepos
    set +x
fi

    
driftedname="$(basename $deposout .npz)"
pvdrifted="$pvdir/${driftedname}.vtp"
speed="1.56*mm/us"
if [ ! -f "$pvdrifted" ] ; then
    echo "FIXME: hard wired speed of $speed"
    set -x
    wirecell-img paraview-depos --speed "$speed" $deposout $pvdrifted
    set +x
fi


echo "Dump blobs to paraview and make blob plots"
for one in $wctdir/clusters-img-*.$clustersext
do
    base="$(basename $one .clusteresext)"
    out="$pvdir/${base}.vtu"
    if [ ! -f "$out" ] ; then
        set -x
        wirecell-img paraview-blobs --speed "$speed" $one $out
        set +x
    fi
    for letter in y z
    do
        pdf="$pltdir/${base}-blobs-t${letter}.pdf"
        if [ -f "$pdf" ] ; then
            echo "already: $pdf"
        else
            set -x
            wirecell-img plot-blobs -p t${letter} $one $pdf
            set +x
        fi
    done
done

