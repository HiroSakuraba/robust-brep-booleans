#!/bin/bash
# S8A shadow validation: run the ENTIRE suite with _LINE_CAST_SHADOW=True.
# Any two-ray vs signed-line disagreement raises LineCastShadowMismatch
# and fails that file. Logs to shadow_run.log.
cd ~/workspace/brep-gates
PY="PYTHONPATH=src:tests $HOME/workspace/brep-booleans/.venv/bin/python"
LOG=tools/review_probes/shadow_run.log
: > "$LOG"
pass=0; fail=0; failed_files=()
run_one() {
  local f="$1"
  echo "=== $f ===" >> "$LOG"
  if timeout 1200 env $PY -c "
import brepkernel.assembly as _A
_A._LINE_CAST_SHADOW = True
import runpy, sys
sys.argv = ['$f']
runpy.run_path('$f', run_name='__main__')
" >> "$LOG" 2>&1; then
    echo "PASS $f"; pass=$((pass+1))
  else
    echo "FAIL $f"; fail=$((fail+1)); failed_files+=("$f")
  fi
}
for f in tests/test_*.py; do
  case "$f" in
    tests/test_s3_fast_paths.py|tests/test_s4_query_context.py|tests/test_s6_session.py|tests/test_s7_many.py)
      echo "=== $f (pytest) ===" >> "$LOG"
      if timeout 1200 env $PY -c "
import brepkernel.assembly as _A
_A._LINE_CAST_SHADOW = True
import pytest, sys
sys.exit(pytest.main(['-x', '-q', '$f']))
" >> "$LOG" 2>&1; then
        echo "PASS $f"; pass=$((pass+1))
      else
        echo "FAIL $f"; fail=$((fail+1)); failed_files+=("$f")
      fi
      ;;
    *) run_one "$f" ;;
  esac
done
echo "RESULT pass=$pass fail=$fail" | tee -a "$LOG"
printf 'failed: %s\n' "${failed_files[@]}" | tee -a "$LOG"
