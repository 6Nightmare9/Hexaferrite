"""Synthetic pw.x output for tests (only the lines hexw.dft.parse reads)."""


def synthetic_pw_out(energy_Ry: float, moments, *, converged=True, iterations=11, job_done=True, wall="12.00s",
                     total_mag=0.0, ts_Ry=0.0, ns=None) -> str:
    """Minimal pw.x-like output with exactly the lines hexw.dft.parse reads. `ns`: {atom (1-based): (up, down)} -> Tr[ns] lines
    of the final SCF step (QE prints them for the atoms with a Hubbard U)."""
    lines = ["     Program PWSCF v.7.6 starts on 26Oct2026 at  6:20:47 ",
             "     number of atoms/cell      =           %d" % len(moments),
             "     number of electrons       =       100.00",
             "     number of k points=     8",
             "     Estimated max dynamical RAM per process >     203.39 MB"]
    for a, (up, dn) in sorted((ns or {}).items()):
        lines.append("     Tr[ns(%3d)] (up, down, total) =   %.5f  %.5f  %.5f" % (a, up, dn, up + dn))
    lines.append("!    total energy              =  %.8f Ry" % energy_Ry)
    lines.append("     estimated scf accuracy    <       0.00000006 Ry")
    lines.append("     smearing contrib. (-TS)   =   %12.8f Ry" % ts_Ry)
    lines.append("     total magnetization       =  %8.2f Bohr mag/cell" % total_mag)
    lines.append("     absolute magnetization    =    17.00 Bohr mag/cell")
    lines.append("     Magnetic moment per site  (integrated on atomic sphere of radius R)")
    for i, m in enumerate(moments, start=1):
        lines.append("     atom %3d (R=0.148)  charge= 12.9006  magn= %8.4f" % (i, m))
    if converged:
        lines.append("     convergence has been achieved in %2d iterations" % iterations)
    else:
        lines.append("     convergence NOT achieved after 200 iterations: stopping")
    lines.append("     PWSCF        :   %s CPU   %s WALL" % (wall, wall))
    if job_done:
        lines.append("   JOB DONE.")
    return "\n".join(lines) + "\n"
