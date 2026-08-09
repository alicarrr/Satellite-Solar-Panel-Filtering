#!/bin/bash
set -eu -o pipefail

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PROJECT_ROOT="$SCRIPTDIR/../"
EXDIR="${PROJECT_ROOT}/examples"

# Global state vars
ERRORMESSAGE=""
failure=0

function show_help {
  echo "script to run the libRadtran optical properties tests"
  echo "  Usage: run_optical_properties_test.sh [OPTION]"
  echo ""
  echo "  -t, --test        <testname>    run a specific test which matches the given substring"
  echo "  -g, --gen_target                (re-)generate optical properties target against which we check"
  echo "  --databasepath    <dirname>     directory which holds the test target output, e.g. test/libRadtran_test_data"
  echo "  -v, --verbose                   Increase verbosity"
  echo "  -h, --help                      Show this help text and exit"
}

function a_gt_b() {
  awk -v n1="$1" -v n2="$2" 'BEGIN {printf (n1>n2?"1":"0") }'
}

function call_uvspec_and_compare_results
{
  inpstr="$1"
  input="$2"
  num_targets="$3"
  all_targets=("${@:4:$num_targets}")
  i=4
  (( i += $num_targets))
  VERBOSE="${@:$i}"
  ((i++))
  QUIET="${@:$i}"
  ((i++))
  GEN_OUTPUT="${@:$i}"

  inpfile=`basename $input`
  error=0

  $QUIET || printf '%.68s %s' "$inpfile.................................................................... "
  $VERBOSE && (echo "")
  $VERBOSE && (echo "uvspec run test $input, $num_targets targets" >> error_optprop.out)
  cat $input > tmp_optprop.out

  # expected output files
  if (( "${#targets[@]}" == 1 )); then
    all_OUTS=("optical_properties.nc")
  else
    declare -a all_OUTS
    for t in "${all_targets[@]}"
    do
      tb=$(basename $t)
      tb="${tb##*iv}"
      idx="${tb%%.*}"
      all_OUTS+=("optical_properties.iv${idx}.nc")
    done
  fi

  # run uvspec
  BINARY=${PROJECT_ROOT}/bin/uvspec
  $WITH_VALGRIND && BINARY="valgrind $BINARY"
  $VERBOSE && (echo "Running $BINARY < $inpfile")
  $VERBOSE && (echo "$inpstr")
  $VERBOSE && (echo "")
  echo "$inpstr" | ${BINARY} &> tmp.out
  RET=$?
  num_outputs=${#all_OUTS[@]}
  for (( i=0; i<${num_outputs}; i++ ));  
  do
    OUT=${all_OUTS[i]}
    target=${all_targets[i]}
    if [ ! -e "$OUT" ] ; then
      echo "I expected that uvspec would create output file $OUT but it didnt :("
      failure=$(( failure + 1 ))
      cat tmp.out
      ERRORMESSAGE+="ERROR: Test $input failed.\n"
      ERRORMESSAGE+=`cat tmp.out`
      return
    fi

    if [ $RET -ne 0 ]
    then
      failure=$(( failure + 1 ))
      cat tmp.out
      ERRORMESSAGE+="ERROR: Test $input failed.\n"
      ERRORMESSAGE+=`cat tmp.out`
      return
    else
      $GEN_OUTPUT && (cp "$OUT" "$target")
    fi
    cat tmp.out 2>&1 >> tmp_optprop

    $VERBOSE && (echo "Difference in optical properties" >> tmp_optprop.out)

    # Filter ncdump output to get a list of variables only
    VARS=$(ncdump -h "$target" | \
      sed -n '/^variables:/,/^$/p;/^$/q' | \
      awk '{print $2}' | \
      awk -F "(" '{print $1}' | \
      sed 's/\\//g')

    if [ "no" == "yes" ]; then
      $VERBOSE && (echo "Comparing $OUT to $target with ")
      # compare meta data first
      set +e
       -mgfbN --tolerance=$ATOL "$OUT" $target >> tmp_optprop.out 2>&1
      let error+=$?
      set -e

      # then look for diffs in each variable
      for varname in $VARS; do
        set +e
         -v ${varname} -n 2 -d --tolerance=$ATOL "$OUT" $target >> tmp_optprop.out 2>&1
        ATOLERR=$?
        if [ $ATOLERR -ne 0 ]; then
           -v ${varname} -n 2 -dfbN --Tolerance=$RTOL "$OUT" $target 2>&1 >> tmp_optprop.out
          RTOLERR=$?
          if [ $RTOLERR -eq 0 ]; then
            $VERBOSE && (echo "Found difference but it is smaller than the relative tolerance")
          fi
          let error+=$RTOLERR
        fi
        set -e
      done

    elif [ "no" == "yes" ]; then
      $VERBOSE && (echo "Comparing $OUT to $target with ")

      for varname in $VARS
      do
        set +e
        NCDIFF="tmp.ncdiff.$varname"
        NCDIFFVALS="tmp.ncdiff_reductions.$varname"
         --no_tmp_fl -v ${varname} "$OUT" $target -O -o $NCDIFF 2>&1 >> tmp_optprop.out
        retcode=$?
        if [ $retcode -ne 0 ]; then cat tmp_optprop.out; fi
        let error+=$retcode
         -v --no_tmp_fl -s "'${varname}_difference'='${varname}'" $NCDIFF -O -o $NCDIFF;
         -v ${varname} --no_tmp_fl -A $OUT $NCDIFF;
         -v --no_tmp_fl -s "\
          '${varname}_max_diff'='${varname}_difference'.fabs().max();                                   \
          '${varname}_rel_diff'=('${varname}_difference'/(1e-300+'${varname}')).fabs().max()*100;       \
          " $NCDIFF -O -o $NCDIFFVALS 2>&1 >> tmp_optprop.out
        ATOLERR=$(ncdump $NCDIFFVALS | grep "${varname}_max_diff = " | awk '{print $3}')
        if (( $(a_gt_b "$ATOLERR" "$ATOL") )); then
          $VERBOSE && (echo "Found atol differences in variable ${varname} with max difference $ATOLERR")
          # check for relative error
          RTOLERR=$(ncdump $NCDIFFVALS | grep "${varname}_rel_diff = " | awk '{print $3}')
          if (( $(a_gt_b "$RTOLERR" "$RTOL") )); then
            error=$(( error + 1 ))
          else
            $VERBOSE && (echo "Found difference but it is smaller than the relative tolerance")
          fi
        fi
        rm -f $NCDIFF
        rm -f $NCDIFFVALS
        set -e
      done
    fi

    $VERBOSE && (echo "Number of differences: $error")
  done
  if [ $error -ne 0 ]
  then
    failure=$(( failure + 1 ))
    ( echo -e "\t\t ... failed." && cat tmp_optprop.out ) 2>&1 >> error_optprop.out
    ERRORMESSAGE+="ERROR: Test $input failed.\n"
    $QUIET || (echo " failed.")
  else
    echo -e "\t\t ... ok." 2>&1 >> error_optprop.out
    $QUIET || (echo " ok.")
  fi
  rm $OUT
  rm error_optprop.out
  rm tmp_optprop.out
  rm tmp_optprop
  rm tmp.out
}


# Check Arguments
ATOL=1e-6 # absolute tolerance that we aim for
RTOL=1e-4 # relative tolerance for variables [%]
GEN_OUTPUT=false
QUIET=false
VERBOSE=false
TEST=""
DATABASEPATH="${PROJECT_ROOT}/examples/write_optical_properties/"
WITH_VALGRIND=false

POSITIONAL=()
while [[ $# -gt 0 ]]
do
  key="$1"

  case $key in
    --atol)
      ATOL="$2"
      shift # past argument
      shift # past value
      ;;
    --rtol)
      RTOL="$2"
      shift # past argument
      shift # past value
      ;;
    -g|--generate-output)
      GEN_OUTPUT=true
      shift # past argument
      ;;
    -q|--quiet)
      VERBOSE2=true
      shift # past argument
      ;;
    -v|--verbose)
      VERBOSE=true
      shift # past argument
      ;;
    --databasepath)
      DATABASEPATH="$2"
      shift # past argument
      shift # past value
      ;;
    --valgrind)
      WITH_VALGRIND=true
      shift # past argument
      ;;
    -t|--test)
      TEST="$2"
      shift # past argument
      shift # past value
      ;;
    -h|--help)
      shift # past argument
      show_help
      exit 0
      ;;
    *)    # unknown option -- ignore
      shift
      ;;
  esac
done

if $VERBOSE ; then
  echo "Found arguments:"
  echo "--atol $ATOL"
  echo "--rtol $RTOL"
  echo "--test $TEST"
  echo "--verbose $VERBOSE"
  echo "--valgrinad $WITH_VALGRIND"
  echo "--quiet $QUIET"
  echo "--generate-target $GEN_OUTPUT"
  echo "--databasepath $DATABASEPATH"
  echo ""
fi

$VERBOSE && echo "Uvspec run optical properties tests"

$GEN_OUTPUT && ($VERBOSE && echo "Making sure that DATABASEPATH exists... "; mkdir -p $DATABASEPATH)

cd $SCRIPTDIR

for input in ${EXDIR}/write_optical_properties/*.INP;
do
  FILENAME=$(basename ${input%.INP})

  if [[ "${FILENAME}" == *"${TEST}"* ]]
  then
    $VERBOSE && echo "Test ${FILENAME} matches given substring $TEST"
  else
    $VERBOSE && echo "Test ${FILENAME} does not match given substring $TEST"
    continue
  fi

  if [[ "${FILENAME}" == *"OPAC"* ]]; then if [[ "no" != "yes" ]]; then
    $VERBOSE && echo "Skipping Test ${FILENAME} because we do not have OPAC";
    continue;
  fi; fi

  [ ! -e "$input" ] && (echo "input filename $input does not exist!"; exit 1)
  targets=($(ls ${DATABASEPATH}/${FILENAME%.INP}*.OUT.nc))
  (( "${#targets[@]}" == 0 )) && ($GEN_OUTPUT || (echo "no benchmark files do exist!"; exit 1))

  inpstr=`cat $input`

  call_uvspec_and_compare_results "$inpstr" "$input" "${#targets[@]}" "${targets[@]}" "$VERBOSE" "$QUIET" "$GEN_OUTPUT" "$WITH_VALGRIND"
done


$QUIET || echo "$ERRORMESSAGE"
$QUIET || echo "${failure} optprop tests failed"
exit $failure
