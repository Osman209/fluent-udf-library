"""
check_conditional_macros.py

Fluent scans a .c file for DEFINE_ macros to build udf_names before the
compiler runs, so it does not see preprocessor conditions. A DEFINE_
macro inside #if ... #endif therefore gets its name registered and then
fails to link:

    lld-link: error: undefined symbol: z_mom_damp

This catches that before anyone hits it in Fluent. Guard the body of the
function, never the macro.
"""
import re, sys, glob, os

bad = []
root = os.path.join(os.path.dirname(__file__), "..")
for path in sorted(glob.glob(os.path.join(root, "*", "*.c"))):
    depth = 0
    for n, line in enumerate(open(path, encoding="utf-8"), 1):
        s = line.strip()
        if re.match(r"#\s*if", s):
            depth += 1
        elif re.match(r"#\s*endif", s):
            depth = max(0, depth - 1)
        elif depth > 0 and re.match(r"DEFINE_[A-Z_]+\s*\(", s):
            bad.append((os.path.relpath(path, root), n, s[:60]))

if bad:
    print("FAIL: DEFINE_ macro inside a preprocessor condition")
    for f, n, t in bad:
        print(f"  {f}:{n}  {t}")
    sys.exit(1)
print(f"ok: no DEFINE_ macro is hidden behind #if in any source file")
