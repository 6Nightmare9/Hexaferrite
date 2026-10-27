"""Run pw_big.exe with OMP threads, poll peak working set (Windows, ctypes). Usage:
run_timed.py <workdir> <input> <omp_threads> <timeout_s>"""
import ctypes, ctypes.wintypes as w, subprocess, sys, time, os, json
wd, inp, nth, tmo = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
obt = sys.argv[5] if len(sys.argv) > 5 else "1"
QE = r"D:\Work\Hexx\tools\qe_ucrt64\ucrt64\bin"
env = dict(os.environ); env["PATH"] = QE + ";" + env["PATH"]
env["OMP_NUM_THREADS"] = nth; env["OPENBLAS_NUM_THREADS"] = obt; env["OMP_STACKSIZE"] = "256M"
env["ESPRESSO_PSEUDO"] = r"D:\Work\Hexx\tools\pseudo"
class PMC(ctypes.Structure):
    _fields_ = [("cb", w.DWORD), ("PageFaultCount", w.DWORD), ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t), ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPagedPoolUsage", ctypes.c_size_t), ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t), ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t)]
k = ctypes.WinDLL("kernel32", use_last_error=True); ps = ctypes.WinDLL("psapi")
k.OpenProcess.restype = w.HANDLE
out = os.path.splitext(inp)[0] + ".out"
t0 = time.time(); peak = 0; status = "done"
with open(os.path.join(wd, inp)) as fi, open(os.path.join(wd, out), "w") as fo:
    p = subprocess.Popen([os.path.join(QE, "pw_big.exe")], stdin=fi, stdout=fo, stderr=subprocess.STDOUT, cwd=wd, env=env)
    h = k.OpenProcess(0x1000 | 0x0400, False, p.pid)
    m = PMC(); m.cb = ctypes.sizeof(m)
    while p.poll() is None:
        if ps.GetProcessMemoryInfo(h, ctypes.byref(m), m.cb): peak = max(peak, m.PeakWorkingSetSize)
        if time.time() - t0 > tmo:
            p.kill(); status = "timeout"; break
        time.sleep(0.5)
    p.wait()
el = time.time() - t0
r = dict(input=inp, omp_threads=nth, blas_threads=obt, wall_s=round(el, 1), peak_ws_GB=round(peak / 2**30, 2), status=status, rc=p.returncode)
json.dump(r, open(os.path.join(wd, os.path.splitext(inp)[0] + ".timing.json"), "w")); print(r)
