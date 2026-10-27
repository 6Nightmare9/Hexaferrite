# install_log

Все установки только внутри D:\Work\Hexx\tools\, без admin; WSL/Docker/драйверы не трогались.

## Движок
- Quantum ESPRESSO 7.6 (pw.x), пакет MSYS2 UCRT64 `mingw-w64-ucrt-x86_64-quantum-espresso-7.6-1`,
  распакован вручную (.pkg.tar.zst) в `tools\qe_ucrt64\ucrt64\bin` вместе с замыканием зависимостей
  (openblas 0.3.34, fftw 3.3.11, hdf5 2.2.0, libgfortran/libgomp/libgcc 16.2.0-6, zlib-ng 2.3.3, libmbd 0.15.0 и др.);
  список с sha256: `tools\_dl\closure_short.tsv`, скачанные пакеты: `tools\_dl\pkgs`.
- Сборка: "Serial multi-threaded version" (OpenMP есть, MPI НЕТ: MPI-сборки QE для win-64 в MSYS2 нет).
- sha256 бинарников: `tools\qe_ucrt64\binary_sha256.txt`.
- `pw_big.exe` = копия `pw.exe` с увеличенным stack reserve в PE-заголовке (`tools\patch_stack.py`);
  нужен для больших ячеек. Запуск: `dft_feasibility\pilots\run_timed.py` (wall и пиковый working set).
- Кандидаты conda-forge win-64 (gpaw, cp2k, siesta, abinit) и pip gpaw/pyscf: по метаданным предыдущего исполнителя
  пригодных готовых сборок для нужной задачи нет; в этом запуске заново не перепроверялось.

## Проверка работоспособности (этот запуск)
- Пилоты a, b, c отработали pw_big.exe с rc=0. Висящих pw.exe/pw_big.exe от прошлого запуска не было (tasklist);
  после работы тоже нет.
- OpenMP: OMP_NUM_THREADS>=2 в этой сборке замедляет/зависает (env_test*.log). OPENBLAS_NUM_THREADS 4-8
  выигрыша не даёт (~82 с/ит против 79 с на 64 атомах). Итог: 1 поток на процесс, параллелизм только запуском
  независимых процессов.

## Псевдопотенциалы (`tools\pseudo\`)
- Источник: SSSP 1.3.0 PBE efficiency (Materials Cloud), архив `sssp_1.3.0_pbe_eff_dl\SSSP_1.3.0_PBE_efficiency.tar.gz`,
  sha256 7a85b71fa3d68df1b5ed33c55c7057681fbf10377ae3b3eac9392924d9189f12; json 2c8642a70308ecd072a209a30c160bf51423560bb202116cb31fd5275eff814a.
  Лицензия набора CC-BY-4.0; семейство GBRV внутри -- GPL-3.0 (Garrity et al., Comput. Mater. Sci. 81, 446 (2014));
  см. `LICENSE.txt` рядом.
- Используются (`pseudo\use`), sha256:
  - fe_pbe_v1.5.uspp.F.UPF 6c132e6e7ac1c42df5cb6a13d0c99bb7af5bf3d6e2349f5658f0c22abfe3911a
  - o_pbe_v1.2.uspp.F.UPF 0640f61fa934269f2500e1cb91f01e48c1d0c3c4ce924f62c2ab35567cfdea41
  - sr_pbe_v1.uspp.F.UPF 292f7f1da5cae4b18da817ab2e82943ad90e74f267a06a6172db8dbcb6edf0ce
  - Zn_pbe_v1.uspp.F.UPF 7710c776159022cdfa14188f04942712811d9a24ad0e74db5ce7c3bdfb95dbbf
  - Co_pbe_v1.2.uspp.F.UPF a5f392a0d35fc78cc9e2ac1636bb1cd02e9b80a063ea368e85be83b0e4523c45
  - mn_pbe_v1.5.uspp.F.UPF c5b6117e22f571895b69dbec8491df6afb2efee1206dad54680881dbd33eb8f0
  - запасные PAW (pslibrary): Fe.pbe-spn-kjpaw_psl.0.2.1.UPF 2317a0d3b136c09b6263e1ba328b505fb1ebf33bb4b28dc3df9e39af3184ef19,
    O.pbe-n-kjpaw_psl.0.1.UPF 7c4b6ed541f83d0afdf5c1d3a8c611340f073f3e2b11e95b7963bb2ee26929aa
- GBRV рассчитаны на ecutwfc 40 / ecutrho 200-240 Ry.
