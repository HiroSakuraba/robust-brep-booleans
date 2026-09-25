"""G2.2: unit tests for the patch keep table (24 cells).

One test per cell of the keep table specified in the work plan (section
2.2). The table is the ONLY place the keep rules live; these tests pin
each cell so a future edit cannot silently change a rule.

Table (from the plan; keep=True means the patch survives in the result):

  operation    | A: IN, OUT, ON_SAME, ON_OPP | B: IN, OUT, ON_SAME, ON_OPP
  union        | drop, keep, keep, drop      | drop, keep, drop, drop
  intersection | keep, drop, keep, drop      | keep, drop, drop, drop
  difference   | drop, keep, drop, keep      | keep, drop, drop, drop
               (B/IN is kept reversed for difference)

Written before the table implementation (I4); they fail until it lands.
"""
import sys

sys.path.insert(0, "src")

from brepkernel.assembly import keep_patch


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")
    return bool(cond)


# (op, operand, state, expected_keep)
CELLS = [
    # union
    ("union", "A", "IN", False),
    ("union", "A", "OUT", True),
    ("union", "A", "ON_SAME", True),
    ("union", "A", "ON_OPP", False),
    ("union", "B", "IN", False),
    ("union", "B", "OUT", True),
    ("union", "B", "ON_SAME", False),
    ("union", "B", "ON_OPP", False),
    # intersection
    ("intersection", "A", "IN", True),
    ("intersection", "A", "OUT", False),
    ("intersection", "A", "ON_SAME", True),
    ("intersection", "A", "ON_OPP", False),
    ("intersection", "B", "IN", True),
    ("intersection", "B", "OUT", False),
    ("intersection", "B", "ON_SAME", False),
    ("intersection", "B", "ON_OPP", False),
    # difference (A - B)
    ("difference", "A", "IN", False),
    ("difference", "A", "OUT", True),
    ("difference", "A", "ON_SAME", False),
    ("difference", "A", "ON_OPP", True),
    ("difference", "B", "IN", True),
    ("difference", "B", "OUT", False),
    ("difference", "B", "ON_SAME", False),
    ("difference", "B", "ON_OPP", False),
]


def main():
    assert len(CELLS) == 24, f"expected 24 cells, got {len(CELLS)}"
    ok = True
    for op, operand, state, expected in CELLS:
        got = keep_patch(op, operand, state)
        ok &= check(
            f"g2 keep {op}/{operand}/{state}",
            got is expected,
            f"keep_patch -> {got}, expected {expected}")
    # The table must reject unknown inputs loudly, never guess.
    for bad in [("union", "A", "BOGUS"), ("xor", "A", "IN"),
                ("union", "C", "IN")]:
        try:
            keep_patch(*bad)
            ok &= check(f"g2 keep rejects {bad}", False, "no error raised")
        except (ValueError, KeyError):
            ok &= check(f"g2 keep rejects {bad}", True)
    print("\nALL PASS" if ok else "\nSOME FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
