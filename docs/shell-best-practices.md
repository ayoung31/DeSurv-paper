# Shell Script Best Practices

Patterns and pitfalls for shell scripts in this project (Slurm submissions, preflight checks, watchdog, etc.).

## Directory Handling (NEVER hardcode paths)

```bash
# WRONG - breaks on other machines
cd /home/naimrashid/Downloads/DeSurv-paper

# RIGHT - portable
cd "${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}"
```

## External Command Guards (prevent `set -e` crashes)

```bash
# WRONG - crashes if sacct unavailable
failed=$(sacct -u $USER ...)

# RIGHT - graceful fallback
if command -v sacct &> /dev/null; then
    failed=$(sacct -u $USER ... || true)
else
    echo "sacct not available"
fi
```

## Multi-line Command Output (prevent arithmetic errors)

```bash
# WRONG - sinfo may return multiple lines
CPUS=$(sinfo -h -o "%C")

# RIGHT - take first line
CPUS=$(sinfo -h -o "%C" | head -1)
```

## Process Elapsed Time (use etime, not CPU time)

```bash
# WRONG - $10 in ps aux is CPU time
ps aux | awk '{print $10}'

# RIGHT - etime shows wall-clock elapsed
ps -eo pid,etime,args
```

## Fallible Commands with `set -e` (grep exits 1 on no match)

```bash
# WRONG - script exits if no matches
set -e
matches=$(grep "pattern" file)

# RIGHT - suppress exit code with || true
matches=$(grep "pattern" file || true)

# Also RIGHT - use if statement
if grep -q "pattern" file; then
    # handle match
fi
```

## Environment Differences (local vs HPC)

```bash
# WRONG - assumes module system exists
module load r/4.4.0

# RIGHT - graceful fallback
module load r/4.4.0 2>/dev/null || true

# For critical dependencies, check explicitly
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found in PATH"
    exit 1
fi
```

## Parsing Multi-value Output (avoid word splitting issues)

```bash
# WRONG - breaks on spaces in values
for val in $(sinfo -h -o "%P %a"); do ...

# RIGHT - use read with IFS
sinfo -h -o "%P|%a" | while IFS='|' read partition avail; do
    echo "Partition: $partition, Available: $avail"
done
```

## `local` Keyword Only Works Inside Functions

```bash
# WRONG - 'local' at top level or in while loop causes error with set -e
while read line; do
    local var=0    # ERROR: local only valid in function
done

# RIGHT - use plain variables outside functions
while read line; do
    var=0          # OK at any scope
done
```

## Capturing Exit Status with `|| true` (masks failures)

```bash
# WRONG - || true makes $? always 0
output=$(some_command || true)
if [ $? -ne 0 ]; then  # Never triggers!
    echo "Failed"
fi

# RIGHT - capture status before || true
output=$(some_command 2>&1)
status=$?
if [ $status -ne 0 ]; then
    echo "Failed with exit $status"
fi
# Now safe to continue even if failed
```
