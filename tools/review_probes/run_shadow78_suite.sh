#!/bin/bash
# S8A shadow validation on OCCT 7.8: classifier-heavy subset with
# _LINE_CAST_SHADOW=True. Any two-ray vs signed-line disagreement
# raises LineCastShadowMismatch and fails that file.
cd ~/workspace/brep-gates
PY="PYTHONPATH=src:tests $HOME/workspace/brep-ocp78-venv/bin/python"
LOG=tools/review_probes/shadow78_run.log
: > "$LOG"
pass=0; fail=0; failed_files=()
for f in tests/test_g5_classifier_independence.py \
         tests/test_sx_point_verdict.py \
         tests/test_s4_query_context.py \
         tests/test_brep_pipeline.py \
         tests/test_c2_independent_frames.py \
         tests/test_regression.py \
         tests/test_degenerate.py \
         tests/test_g7_seam_stress.py; do
  echo "=== $f ===" >> "$LOG"
  if timeout 1800 env $PY -c "
import brepkernel.assembly as _A
_A._LINE_CAST_SHADOW = True
import runpy
runpy.run_path('$f', run_name='__main__')
" >> "$LOG" 2>&1; then
    echo "PASS $f"; pass=$((pass+1))
  else
    echo "FAIL $f"; fail=$((fail+1)); failed_files+=("$f")
  fi
done
echo "RESULT pass=$pass fail=$fail" | tee -a "$LOG"
printf 'failed: %s\n' "${failed_files[@]}" | tee -a "$LOG"
