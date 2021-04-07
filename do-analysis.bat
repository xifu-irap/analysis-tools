echo Set Python environment variables...

@echo off
call C:\GNU\WPy-3720\scripts\env.bat

"C:\GNU\WPy64-3720\python-3.7.2.amd64\python.exe"  analyse_demux_perfs.py
"C:\GNU\WPy64-3720\python-3.7.2.amd64\python.exe"  toDACtest.py
