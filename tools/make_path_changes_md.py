import csv, collections, re, pathlib
rows=list(csv.DictReader(open('tools/path_changes.tsv',encoding='utf-8'),delimiter='\t',quoting=csv.QUOTE_NONE))
pat=collections.Counter(); perfile=collections.defaultdict(list)
for r in rows:
    f=r['file']; old=r['old']; new=r['new']; ext=pathlib.Path(f).suffix
    m=re.match(r'^([rRfF]*)(["\'])(G:[\\/]+Comorbidity|/workspace/Comorbidity)([\\/]+([A-Za-z_]+))?',old)
    pre,q,root,_,sub=m.groups()
    key=(ext, f'{pre}{q}{root}' + (('\\' if '\\' in root else '/') + f'{sub}…' if sub else '') + q,
         'parallel worker (clusterEvalQ)' if 'Sys.getenv' in new else ('relocated helper' if 'helpers' in new else 'main process'))
    pat[key]+=1; perfile[f].append(r)
L=[]
L.append("# PATH_CHANGES — hard-coded paths replaced in the copied scripts\n\n")
L.append("Only the copies in `code_release/scripts/` were edited; the originals in the project folder are unchanged. "
         "The edits are mechanical (string-literal substitution plus a short header); no analysis logic was changed. "
         "They were made by `tools/make_release.py`; the line-by-line log (file, line number in the original, old literal, new expression) is `tools/path_changes.tsv`.\n\n")
L.append(f"**Total: {len(rows)} replacements in {len(perfile)} files.**\n\n")
L.append("## 1. Header added to every edited file\n\n")
L.append("R (inserted after the leading comment block):\n\n```r\n# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder\nPROJ_ROOT <- local({\n  a <- sub(\"^--file=\", \"\", grep(\"^--file=\", commandArgs(FALSE), value = TRUE))\n  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()\n  while (basename(d) != \"scripts\" && dirname(d) != d) d <- dirname(d)\n  gsub(\"\\\\\\\\\", \"/\", Sys.getenv(\"COMORBIDITY_ROOT\", unset = normalizePath(if (basename(d) == \"scripts\") dirname(d) else \"..\", winslash = \"/\", mustWork = FALSE)))\n})\nSys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers\n```\n\n")
L.append("Python (inserted after the module docstring):\n\n```python\nimport os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder\nPROJ_ROOT = (_os.environ.get(\"COMORBIDITY_ROOT\") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == \"scripts\"), str(_pl.Path.cwd()))).replace(\"\\\\\", \"/\")\n```\n\n")
L.append("`PROJ_ROOT` uses forward slashes on all platforms. With `Rscript scripts/xx.R` or `python scripts/xx.py` the default resolves to the repository root (the parent of `scripts/`), wherever the command is started from. Both headers were tested (R 4.5.2 on Windows; Python 3.12 on Linux), with and without `COMORBIDITY_ROOT` set.\n\n")
L.append("## 2. Substitution rules\n\n")
L.append("| Language | Original literal | Replacement |\n|---|---|---|\n"
 "| R | `\"G:/Comorbidity\"` | `PROJ_ROOT` |\n"
 "| R | `\"G:/Comorbidity/<sub/path>\"` | `file.path(PROJ_ROOT, \"<sub/path>\")` |\n"
 "| R, inside `clusterEvalQ(...)` (evaluated on parallel workers, where `PROJ_ROOT` is not defined) | `\"G:/Comorbidity/Rlib\"` | `file.path(Sys.getenv(\"COMORBIDITY_ROOT\"), \"Rlib\")` |\n"
 "| R (helper moved into the release) | `\"G:/Comorbidity/tmp/revmr/cause_fix.R\"` | `file.path(PROJ_ROOT, \"scripts/helpers/cause_fix.R\")` |\n"
 "| Python | `r\"G:\\Comorbidity\"`, `\"G:/Comorbidity\"` | `PROJ_ROOT` (e.g. `Path(r\"G:\\Comorbidity\")` → `Path(PROJ_ROOT)`) |\n"
 "| Python | `r\"G:\\Comorbidity\\<sub\\path>\"`, `\"G:/Comorbidity/<sub/path>\"`, `'G:/Comorbidity/<sub/path>'` | `f\"{PROJ_ROOT}/<sub/path>\"` (backslashes converted to `/`) |\n"
 "| Python | `'/workspace/Comorbidity'` (Linux sandbox mount; `revision3/r3_eye_build.py`) | `PROJ_ROOT` |\n\n")
L.append("## 3. Pattern counts\n\n")
L.append("| File type | Original pattern (root + first sub-folder) | Context | n |\n|---|---|---|---|\n")
for (ext,p,ctx),n in sorted(pat.items(), key=lambda x:(x[0][0],-x[1])):
    L.append(f"| {ext} | `{p}` | {ctx} | {n} |\n")
L.append("\n## 4. Replacements per file\n\n")
L.append("| File (under code_release/) | n | Lines in original |\n|---|---|---|\n")
for f in sorted(perfile):
    ls=sorted({int(r['line_in_original']) for r in perfile[f]})
    L.append(f"| `{f}` | {len(perfile[f])} | {', '.join(map(str,ls))} |\n")
allf=sorted(p.as_posix() for p in pathlib.Path('scripts').rglob('*') if p.suffix in ('.R','.py'))
unchanged=[f for f in allf if f not in perfile]
L.append("\n## 5. Copied scripts without hard-coded paths (not edited)\n\n")
L.append("These already use relative paths (`Path(__file__).resolve().parents[1]` in `01_download_nhanes.py` and `02_nhanes_codebook.py`; `Path(__file__).parent` in `build_docx.py` and `build_af1.py`, which expect the manuscript markdown parts next to them), take paths as command-line arguments, or import their paths from `rev_common.py`:\n\n")
for f in unchanged: L.append(f"- `{f}`\n")
L.append("\n## 6. Deliberately not changed\n\n"
 "- Executable names. PLINK, SMR and GCTA are called as Windows binaries under `<root>/data/raw/tools/` (`plink/plink.exe`, `smr-1.3.1-win-x86_64/smr-1.3.1-win.exe`, `gcta*.exe`). On Linux or macOS, put the corresponding binaries there under the same names or edit the `PLINK`/`SMR`/`GCTA` variable definitions (see `environment/external_tools.md`).\n"
 "- Relative sub-paths below the root (`data/raw/...`, `derived/...`, `results/...`, `tmp/...`, `submission/...`). Outputs are written to the same relative locations as in the original project.\n"
 "- `.libPaths(c(file.path(PROJ_ROOT, \"Rlib\"), ...))`. The project-local R library is still prepended; if `<root>/Rlib` does not exist R ignores it and uses the default libraries.\n")
open('PATH_CHANGES.md','w',encoding='utf-8').write(''.join(L))
print(len(rows), len(perfile), len(unchanged))
