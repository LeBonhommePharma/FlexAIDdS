#!/bin/bash
# Answer "did the workers actually start?" from LOG CONTENT, never from file existence.
#
# WHY THIS EXISTS -- measured 2026-09-21
#   Three consecutive campaign relaunches were reported as running on the strength
#   of (a) a launch-receipt JSON written before dispatch and (b) worker log files
#   existing. Every one had died instantly. The logs were 41 bytes:
#       nohup: setsid: No such file or directory
#   setsid is a Linux utility and does not exist on macOS. The receipt recorded
#   INTENT; the log files recorded a failure; neither existence check could tell
#   the difference. The same shape had already cost this project nine campaign
#   cells (an engine that exited rc=2 with a message that read like a docking
#   result) and a test suite that "passed" after exiting in 1 s on an unknown flag.
#
# CONTRACT: exit 0 only when at least one worker log contains POSITIVE evidence of
# the launcher running. Exit 1 on a known failure signature. Exit 3 when there is
# nothing to judge -- an absent or empty log must never read as success.
set -u
LOGDIR=${1:?usage: verify_launch.sh <logdir> <glob> [min_bytes]}
GLOB=${2:?}
MIN=${3:-200}

# Signatures of a dispatch that never became a process. Each was observed.
FAIL_RE='No such file or directory|command not found|Permission denied|cannot execute'
# Evidence the launcher body actually ran (its own output, not the shell's).
OK_RE='PREFLIGHT|\[[0-9]+/[0-9]+\]|skip .* \(done\)|GATES OK|launch window'

shopt -s nullglob
logs=( "$LOGDIR"/$GLOB )
if [ ${#logs[@]} -eq 0 ]; then
  echo "VACUOUS: no logs matching $GLOB in $LOGDIR -- cannot conclude a launch happened" >&2
  exit 3
fi

alive=0 dead=0 empty=0
for f in "${logs[@]}"; do
  b=$(wc -c < "$f" | tr -d ' ')
  if grep -qE "$FAIL_RE" "$f" 2>/dev/null; then
    printf "  FAILED   %-34s %6s B  %s\n" "$(basename "$f")" "$b" "$(head -1 "$f" | cut -c1-58)"
    dead=$((dead+1))
  elif [ "$b" -lt "$MIN" ] || ! grep -qE "$OK_RE" "$f" 2>/dev/null; then
    printf "  NO-PROOF %-34s %6s B  (no launcher output yet)\n" "$(basename "$f")" "$b"
    empty=$((empty+1))
  else
    printf "  RUNNING  %-34s %6s B  %s\n" "$(basename "$f")" "$b" \
      "$(grep -oE "$OK_RE" "$f" | tail -1)"
    alive=$((alive+1))
  fi
done

echo "  alive=$alive failed=$dead no-proof=$empty of ${#logs[@]}"
[ "$dead" -gt 0 ] && { echo "VERDICT: at least one worker never started" >&2; exit 1; }
[ "$alive" -eq 0 ] && { echo "VERDICT: no worker has produced launcher output -- not launched" >&2; exit 3; }
echo "VERDICT: launched"
exit 0
