-- Lmod modulefile template. Point PREFIX at an installed FlexAIDdS prefix.
--   module use $PWD/packaging/modulefiles
--   export FLEXAIDDS_PREFIX=/path/to/install
help([[
FlexAID∆S native engine + optional Python venv.

The native binary and the flexaidds Python package are separate.
]])

local prefix = os.getenv("FLEXAIDDS_PREFIX") or "/usr/local"
prepend_path("PATH", pathJoin(prefix, "bin"))
setenv("FLEXAIDDS_DATA_DIR", pathJoin(prefix, "share"))

local venv = os.getenv("FLEXAIDDS_VENV") or pathJoin(os.getenv("HOME") or "", ".flexaidds/venv")
if isDir(pathJoin(venv, "bin")) then
  prepend_path("PATH", pathJoin(venv, "bin"))
end

whatis("FlexAID∆S docking engine")
