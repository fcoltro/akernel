#!/usr/bin/env bash
# Run ak_random_models over a range of seeds and survive crashes.
#
#   fuzz/nightly.sh EXE OUTDIR FIRST_SEED COUNT SECONDS [more ak_random_models options]
#
# A crash (a signal, an access violation) ends the process, not the run: the
# seed it was on is read from the progress file, recorded in
# OUTDIR/crashes.txt, and the run goes on from the next seed. Each stretch of
# seeds between crashes writes its own report (OUTDIR/part-N.md, .json, .log);
# OUTDIR/summary.md joins them. Exit 1 when any oracle failed or anything
# crashed, else 0.
set -u
exe=$1 out=$2 first=$3 count=$4 seconds=$5
shift 5
mkdir -p "$out"
: > "$out/crashes.txt"
deadline=$(( $(date +%s) + seconds ))
seed=$first end=$(( first + count )) part=0 status=0
while [ "$seed" -lt "$end" ]; do
  left=$(( deadline - $(date +%s) ))
  [ "$left" -le 0 ] && break
  rm -f "$out/running"
  "$exe" --seed "$seed" --count $(( end - seed )) --time-limit "$left" \
         --json "$out/part-$part.json" --markdown "$out/part-$part.md" --progress "$out/running" "$@" \
         > "$out/part-$part.log" 2>&1
  code=$?
  if [ "$code" -eq 0 ]; then break; fi
  if [ "$code" -eq 1 ]; then status=1; break; fi
  status=1
  if [ ! -s "$out/running" ]; then
    echo "exit code $code with no running seed (see part-$part.log)" | tee -a "$out/crashes.txt"
    break
  fi
  crashed=$(tr -d '\r\n ' < "$out/running")
  echo "seed $crashed: the process ended with exit code $code (replay: ak_random_models --seed $crashed --count 1 --verbose)" \
    | tee -a "$out/crashes.txt"
  tail -n 20 "$out/part-$part.log"
  seed=$(( crashed + 1 ))
  part=$(( part + 1 ))
done
{
  for f in "$out"/part-*.md; do [ -f "$f" ] && { cat "$f"; echo; }; done
  if [ -s "$out/crashes.txt" ]; then
    echo "### Crashes"; echo
    sed 's/^/- /' "$out/crashes.txt"
  fi
} > "$out/summary.md"
cat "$out/summary.md"
exit $status
