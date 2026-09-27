#!/bin/bash
# S8A shadow validation (8.0.1): classifier-heavy + geometry-diverse
# subset with _LINE_CAST_SHADOW=True. Any two-ray vs signed-line
# disagreement raises LineCastShadowMismatch and fails that file.
cd ~/workspace/brep-gates
PY="PYTHONPATH=src:tests $HOME/workspace/brep-booleans/.venv/bin/python"
LOG=tools/review_probes/shadow81_run.log
: > "$LOG"
pass=0; fail=0; failed_files=()
run_one() {
  local f="$1"; local pytest_mode="$2"
  echo "=== $f ===" >> "$LOG"
  if [ "$pytest_mode" = "pytest" ]; then
    CMD="import brepkernel.assembly as _A; _A._LINE_CAST_SHADOW = True
import pytest, sys; sys.exit(pytest.main(['-x', '-q', '$f']))"
  else
    CMD="import brepkernel.assembly as _A; _A._LINE_CAST_SHADOW = True
import runpy; runpy.run_path('$f', run_name='__main__')"
  fi
  if timeout 2400 env $PY -c "$CMD" >> "$LOG" 2>&1; then
    echo "PASS $f"; pass=$((pass+1))
  else
    echo "FAIL $f"; fail=$((fail+1)); failed_files+=("$f")
  fi
}
run_one tests/test_g5_classifier_independence.py
run_one tests/test_sx_point_verdict.py
run_one tests/test_s4_query_context.py pytest
run_one tests/test_brep_pipeline.py
run_one tests/test_c2_independent_frames.py
run_one tests/test_c3_self_touch.py
run_one tests/test_regression.py
run_one tests/test_degenerate.py
run_one tests/test_g7_seam_stress.py
run_one tests/test_s8a_ray_avoidance.py pytest
run_one tests/test_freeform_nurbs.py
run_one tests/test_rigid_motion_invariance.py
echo "RESULT pass=$pass fail=$fail" | tee -a "$LOG"
printf 'failed: %s\n' "${failed_files[@]}" | tee -a "$LOG"
