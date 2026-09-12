"""
check_udf_patterns.py

Static checks for the mistakes that only show up in a real Fluent run,
each one added after it was actually hit:

  1. A DEFINE_ macro inside #if. Fluent scans for DEFINE_ macros to build
     udf_names before the compiler runs, registers the name, and then
     fails to link it.

  2. A "have I done this yet" static flag. A loaded library keeps its
     statics across Initialize, so on a second run the flag never fires:
     logs get appended to instead of rewritten, "log only when the clock
     advanced" tests stop passing, and captured reference positions
     silently belong to the previous run. Use udf_restarted instead.

  3. Writing a file without UDF_IS_WRITER, so every process writes.
     (Reading is fine on every process and is not flagged.)

  4. Accumulating over faces without PRINCIPAL_FACE_P, or over cells with
     begin_c_loop instead of begin_c_loop_int, both of which double count
     in parallel.

Only the UDF packs are checked; tests/ is ordinary C run with gcc.
"""
import re, sys, glob, os

ROOT = os.path.join(os.path.dirname(__file__), "..")
PACKS = ("motion", "wavetank", "profiles")
problems = []

def flag(path, line, msg):
    problems.append(f"{os.path.relpath(path, ROOT)}:{line}  {msg}")

for pack in PACKS:
    for path in sorted(glob.glob(os.path.join(ROOT, pack, "*.c"))):
        lines = open(path, encoding="utf-8").read().split("\n")

        depth = 0
        for i, raw in enumerate(lines, 1):
            s = raw.strip()

            # 1. DEFINE_ behind a preprocessor condition
            if re.match(r"#\s*if", s):
                depth += 1
            elif re.match(r"#\s*endif", s):
                depth = max(0, depth - 1)
            elif depth > 0 and re.match(r"DEFINE_[A-Z_]+\s*\(", s):
                flag(path, i, "DEFINE_ macro inside #if: registers but will not link")

            # 2. first-call flag
            if re.search(r"static\s+(int|short|char)\s+first\b", s):
                flag(path, i, "static first-call flag survives Initialize: use udf_restarted")

            # 3. writing a file from every process
            m = re.search(r'fopen\s*\([^,]+,\s*"([aw])', s)
            if m:
                ctx = "\n".join(lines[max(0, i - 30):i])
                if "UDF_IS_WRITER" not in ctx:
                    flag(path, i, "file opened for writing without UDF_IS_WRITER")

            # 4. parallel double counting
            if re.match(r"begin_f_loop\s*\(", s):
                ctx = "\n".join(lines[i - 1:i + 30])
                if re.search(r"\w+\s*\[[^\]]*\]\s*\+=", ctx) and "PRINCIPAL_FACE_P" not in ctx:
                    flag(path, i, "face loop accumulates without PRINCIPAL_FACE_P")
            if re.match(r"begin_c_loop\s*\(", s):
                ctx = "\n".join(lines[i - 1:i + 30])
                if re.search(r"\w+\s*\[[^\]]*\]\s*\+=", ctx):
                    flag(path, i, "begin_c_loop accumulates: use begin_c_loop_int")

if problems:
    print("FAIL: patterns that break in a real Fluent run")
    for p in problems:
        print("  " + p)
    sys.exit(1)

n = sum(len(glob.glob(os.path.join(ROOT, p, "*.c"))) for p in PACKS)
print(f"ok: {n} UDF files clear of the four patterns that only fail inside Fluent")
