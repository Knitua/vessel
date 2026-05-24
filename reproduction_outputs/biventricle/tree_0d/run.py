import atexit
import os
import signal
import subprocess
import sys

PINNED_SOLVER = '/root/autodl-tmp/svVascularize/svv/utils/solvers/0D/Linux/x86_64/svzerodsolver'
INPUT_FILENAME = 'solver_0d.in'
CHILD_PID_FILE = ".svv_0d_solver.pid"
_CHILD_PROCESS = None


def _pid_path(run_dir):
    return os.path.join(run_dir, CHILD_PID_FILE)


def _write_child_pid(run_dir, pid):
    try:
        with open(_pid_path(run_dir), "w", encoding="utf-8") as pid_file:
            pid_file.write(str(pid))
    except Exception:
        pass


def _remove_child_pid(run_dir):
    try:
        os.remove(_pid_path(run_dir))
    except FileNotFoundError:
        pass
    except Exception:
        pass


def _resolve_solver():
    if PINNED_SOLVER and os.path.isfile(PINNED_SOLVER):
        return PINNED_SOLVER
    try:
        from svv.utils.solvers.solver_0d import get_solver_0d_exe

        solver = str(get_solver_0d_exe())
        if os.path.isfile(solver):
            return solver
    except Exception:
        pass
    return None


def _print_setup_steps():
    print("No svZeroDSolver executable is available.", file=sys.stderr, flush=True)
    print("Expected a packaged or locally built 0D solver from svVascularize.", file=sys.stderr, flush=True)
    print("Setup options:", file=sys.stderr, flush=True)
    print("  1) Build 0D solver binary:", file=sys.stderr, flush=True)
    print("     python setup.py build_ext --build-solver-0d", file=sys.stderr, flush=True)
    print("  2) Build all native binaries (MMG + 0D):", file=sys.stderr, flush=True)
    print("     python setup.py build_ext --build-native-binaries", file=sys.stderr, flush=True)
    print("  3) Use explicit executable path:", file=sys.stderr, flush=True)
    print("     export SVV_SOLVER_0D_PATH=/abs/path/to/svzerodsolver", file=sys.stderr, flush=True)
    print("  4) Re-export with path_to_0d_solver='/abs/path/to/svzerodsolver'", file=sys.stderr, flush=True)


def _stop_child(run_dir, force=False):
    global _CHILD_PROCESS
    proc = _CHILD_PROCESS
    if proc is None:
        _remove_child_pid(run_dir)
        return
    if proc.poll() is not None:
        _remove_child_pid(run_dir)
        return
    try:
        if force:
            proc.kill()
        else:
            proc.terminate()
            try:
                proc.wait(timeout=3.0)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait(timeout=3.0)
    except Exception:
        pass
    finally:
        _remove_child_pid(run_dir)


def _handle_signal(signum, _frame):
    run_dir = os.path.dirname(os.path.abspath(__file__))
    signame = str(signum)
    try:
        signame = signal.Signals(signum).name
    except Exception:
        pass
    print("Received " + signame + "; stopping svZeroDSolver...", file=sys.stderr, flush=True)
    _stop_child(run_dir, force=False)
    raise SystemExit(130)


def main():
    run_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(run_dir, INPUT_FILENAME)
    if not os.path.isfile(input_file):
        print("Input file not found: " + input_file, file=sys.stderr, flush=True)
        sys.exit(2)

    solver = _resolve_solver()
    if solver is None:
        _print_setup_steps()
        sys.exit(2)

    for signame in ("SIGTERM", "SIGINT", "SIGHUP"):
        if hasattr(signal, signame):
            signal.signal(getattr(signal, signame), _handle_signal)

    atexit.register(_stop_child, run_dir, False)

    print("Running svZeroDSolver: " + solver, flush=True)
    global _CHILD_PROCESS
    _CHILD_PROCESS = subprocess.Popen([solver, input_file], cwd=run_dir)
    _write_child_pid(run_dir, _CHILD_PROCESS.pid)
    try:
        return_code = _CHILD_PROCESS.wait()
    finally:
        _remove_child_pid(run_dir)
    sys.exit(return_code)


if __name__ == "__main__":
    main()
